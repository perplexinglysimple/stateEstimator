#include "EKF.h"

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define EARTH_RADIUS_M 6378137.0
#define GRAVITY_MPS2 9.80665
#define PI 3.14159265358979323846
#define DEFAULT_PERIOD_MS 1000LL
#define IMU_FRESHNESS_MS 5000LL
#define GPS_FRESHNESS_MS 30000LL

typedef struct
{
    long long tick_ms;
    double accel_x_mps2;
    double accel_y_mps2;
    double accel_z_mps2;
} ImuRow;

typedef struct
{
    long long tick_ms;
    int valid_fix;
    int fix_type;
    int satellites;
    double latitude_deg;
    double longitude_deg;
    double altitude_m;
    double speed_mps;
    double course_deg;
    double hdop;
} GpsRow;

typedef struct
{
    long long* values;
    size_t len;
    size_t cap;
} TickArray;

typedef struct
{
    ImuRow* values;
    size_t len;
    size_t cap;
} ImuArray;

typedef struct
{
    GpsRow* values;
    size_t len;
    size_t cap;
} GpsArray;

typedef struct
{
    double dt;
    bool have_accel;
    double accel[3];
} PredictContext;

static int append_tick(TickArray* arr, long long value)
{
    long long* newBuf;
    size_t newCap;
    if (arr->len < arr->cap)
    {
        arr->values[arr->len++] = value;
        return 0;
    }
    newCap = arr->cap == 0 ? 1024 : arr->cap * 2;
    newBuf = (long long*) realloc(arr->values, newCap * sizeof(long long));
    if (newBuf == NULL)
    {
        return -1;
    }
    arr->values = newBuf;
    arr->cap = newCap;
    arr->values[arr->len++] = value;
    return 0;
}

static int append_imu(ImuArray* arr, const ImuRow* row)
{
    ImuRow* newBuf;
    size_t newCap;
    if (arr->len < arr->cap)
    {
        arr->values[arr->len++] = *row;
        return 0;
    }
    newCap = arr->cap == 0 ? 1024 : arr->cap * 2;
    newBuf = (ImuRow*) realloc(arr->values, newCap * sizeof(ImuRow));
    if (newBuf == NULL)
    {
        return -1;
    }
    arr->values = newBuf;
    arr->cap = newCap;
    arr->values[arr->len++] = *row;
    return 0;
}

static int append_gps(GpsArray* arr, const GpsRow* row)
{
    GpsRow* newBuf;
    size_t newCap;
    if (arr->len < arr->cap)
    {
        arr->values[arr->len++] = *row;
        return 0;
    }
    newCap = arr->cap == 0 ? 1024 : arr->cap * 2;
    newBuf = (GpsRow*) realloc(arr->values, newCap * sizeof(GpsRow));
    if (newBuf == NULL)
    {
        return -1;
    }
    arr->values = newBuf;
    arr->cap = newCap;
    arr->values[arr->len++] = *row;
    return 0;
}

static int compare_ll(const void* a, const void* b)
{
    const long long va = *(const long long*) a;
    const long long vb = *(const long long*) b;
    if (va < vb)
    {
        return -1;
    }
    if (va > vb)
    {
        return 1;
    }
    return 0;
}

static int compare_imu_by_tick(const void* a, const void* b)
{
    const ImuRow* ra = (const ImuRow*) a;
    const ImuRow* rb = (const ImuRow*) b;
    if (ra->tick_ms < rb->tick_ms)
    {
        return -1;
    }
    if (ra->tick_ms > rb->tick_ms)
    {
        return 1;
    }
    return 0;
}

static int compare_gps_by_tick(const void* a, const void* b)
{
    const GpsRow* ra = (const GpsRow*) a;
    const GpsRow* rb = (const GpsRow*) b;
    if (ra->tick_ms < rb->tick_ms)
    {
        return -1;
    }
    if (ra->tick_ms > rb->tick_ms)
    {
        return 1;
    }
    return 0;
}

static bool parse_long_long(const char* text, long long* out)
{
    char* endptr = NULL;
    long long v;
    if (text == NULL || *text == '\0')
    {
        return false;
    }
    errno = 0;
    v = strtoll(text, &endptr, 10);
    if (errno != 0 || endptr == text)
    {
        return false;
    }
    *out = v;
    return true;
}

static bool parse_int(const char* text, int* out)
{
    char* endptr = NULL;
    long v;
    if (text == NULL || *text == '\0')
    {
        return false;
    }
    errno = 0;
    v = strtol(text, &endptr, 10);
    if (errno != 0 || endptr == text)
    {
        return false;
    }
    *out = (int) v;
    return true;
}

static bool parse_double(const char* text, double* out)
{
    char* endptr = NULL;
    double v;
    if (text == NULL || *text == '\0')
    {
        return false;
    }
    errno = 0;
    v = strtod(text, &endptr);
    if (errno != 0 || endptr == text)
    {
        return false;
    }
    *out = v;
    return true;
}

static int split_csv(char* line, char** fields, int maxFields)
{
    int count = 0;
    char* token;
    char* ctx = NULL;
#ifdef _WIN32
    token = strtok_s(line, ",\r\n", &ctx);
#else
    token = strtok_r(line, ",\r\n", &ctx);
#endif
    while (token != NULL && count < maxFields)
    {
        fields[count++] = token;
#ifdef _WIN32
        token = strtok_s(NULL, ",\r\n", &ctx);
#else
        token = strtok_r(NULL, ",\r\n", &ctx);
#endif
    }
    return count;
}

static void join_path(char* out, size_t outSize, const char* dir, const char* file)
{
    size_t len = strlen(dir);
    if (len > 0 && (dir[len - 1] == '\\' || dir[len - 1] == '/'))
    {
        snprintf(out, outSize, "%s%s", dir, file);
    }
    else
    {
        snprintf(out, outSize, "%s\\%s", dir, file);
    }
}

static FILE* open_read_file(const char* filepath)
{
    FILE* f = NULL;
#ifdef _WIN32
    if (fopen_s(&f, filepath, "r") != 0)
    {
        f = NULL;
    }
#else
    f = fopen(filepath, "r");
#endif
    return f;
}

static FILE* open_write_file(const char* filepath)
{
    FILE* f = NULL;
#ifdef _WIN32
    if (fopen_s(&f, filepath, "w") != 0)
    {
        f = NULL;
    }
#else
    f = fopen(filepath, "w");
#endif
    return f;
}

static int load_imu_rows(const char* filepath, ImuArray* imu, TickArray* ticks)
{
    FILE* f = open_read_file(filepath);
    char line[1024];
    if (f == NULL)
    {
        fprintf(stderr, "Failed to open IMU CSV: %s\n", filepath);
        return -1;
    }
    if (fgets(line, sizeof(line), f) == NULL)
    {
        fclose(f);
        return 0;
    }
    while (fgets(line, sizeof(line), f) != NULL)
    {
        char* fields[20];
        ImuRow row;
        int n = split_csv(line, fields, 20);
        if (n < 4)
        {
            continue;
        }
        if (!parse_long_long(fields[0], &row.tick_ms))
        {
            continue;
        }
        if (!parse_double(fields[1], &row.accel_x_mps2))
        {
            continue;
        }
        if (!parse_double(fields[2], &row.accel_y_mps2))
        {
            continue;
        }
        if (!parse_double(fields[3], &row.accel_z_mps2))
        {
            continue;
        }
        if (append_imu(imu, &row) != 0 || append_tick(ticks, row.tick_ms) != 0)
        {
            fclose(f);
            return -1;
        }
    }
    fclose(f);
    return 0;
}

static int load_gps_rows(const char* filepath, GpsArray* gps, TickArray* ticks)
{
    FILE* f = open_read_file(filepath);
    char line[1024];
    if (f == NULL)
    {
        fprintf(stderr, "Failed to open GPS CSV: %s\n", filepath);
        return -1;
    }
    if (fgets(line, sizeof(line), f) == NULL)
    {
        fclose(f);
        return 0;
    }
    while (fgets(line, sizeof(line), f) != NULL)
    {
        char* fields[16];
        GpsRow row;
        int n = split_csv(line, fields, 16);
        if (n < 10)
        {
            continue;
        }
        if (!parse_long_long(fields[0], &row.tick_ms))
        {
            continue;
        }
        if (!parse_int(fields[1], &row.valid_fix))
        {
            continue;
        }
        if (!parse_int(fields[2], &row.fix_type))
        {
            continue;
        }
        if (!parse_int(fields[3], &row.satellites))
        {
            continue;
        }
        if (!parse_double(fields[4], &row.latitude_deg))
        {
            continue;
        }
        if (!parse_double(fields[5], &row.longitude_deg))
        {
            continue;
        }
        if (!parse_double(fields[6], &row.altitude_m))
        {
            continue;
        }
        if (!parse_double(fields[7], &row.speed_mps))
        {
            continue;
        }
        if (!parse_double(fields[8], &row.course_deg))
        {
            continue;
        }
        if (!parse_double(fields[9], &row.hdop))
        {
            continue;
        }
        if (append_gps(gps, &row) != 0 || append_tick(ticks, row.tick_ms) != 0)
        {
            fclose(f);
            return -1;
        }
    }
    fclose(f);
    return 0;
}

static int build_timeline(const TickArray* sortedTicks, long long periodMs, TickArray* timeline)
{
    size_t i;
    long long prev = LLONG_MIN;
    long long nextTick = 0;
    for (i = 0; i < sortedTicks->len; ++i)
    {
        long long tick = sortedTicks->values[i];
        if (tick == prev)
        {
            continue;
        }
        prev = tick;
        if (timeline->len == 0 || tick >= nextTick)
        {
            if (append_tick(timeline, tick) != 0)
            {
                return -1;
            }
            nextTick = tick + periodMs;
        }
    }
    return 0;
}

static bool is_fresh(long long ekfTick, long long sensorTick, long long windowMs)
{
    long long age = ekfTick - sensorTick;
    return age >= 0 && age <= windowMs;
}

static void geodetic_to_local_m(
    double latDeg,
    double lonDeg,
    double altM,
    double refLatDeg,
    double refLonDeg,
    double refAltM,
    double* xEast,
    double* yNorth,
    double* zUp)
{
    double lat = latDeg * PI / 180.0;
    double lon = lonDeg * PI / 180.0;
    double refLat = refLatDeg * PI / 180.0;
    double refLon = refLonDeg * PI / 180.0;
    double dLat = lat - refLat;
    double dLon = lon - refLon;
    *yNorth = dLat * EARTH_RADIUS_M;
    *xEast = dLon * EARTH_RADIUS_M * cos(refLat);
    *zUp = altM - refAltM;
}

static void speed_course_to_local_velocity(double speedMps, double courseDeg, double* vxEast, double* vyNorth, double* vzUp)
{
    double courseRad = courseDeg * PI / 180.0;
    *vyNorth = speedMps * cos(courseRad);
    *vxEast = speedMps * sin(courseRad);
    *vzUp = 0.0;
}

static double clip_value(double value, double low, double high)
{
    if (value < low)
    {
        return low;
    }
    if (value > high)
    {
        return high;
    }
    return value;
}

static void TransitionFunction(EKFMatrix* x, EKFMatrix* xPredicted, EKFState* ekf, void* userData)
{
    int i;
    int j;
    const PredictContext* ctx = (const PredictContext*) userData;
    matrixType pred[6] = {0};
    matrixType dt = (matrixType) ((ctx != NULL) ? ctx->dt : 0.1);

    for (i = 0; i < 6; ++i)
    {
        matrixType sum = 0;
        for (j = 0; j < 6; ++j)
        {
            sum += ACCESS_MATRIX(*(ekf->A), i, j) * ACCESS_MATRIX(*x, j, 0);
        }
        pred[i] = sum;
    }

    if (ctx != NULL && ctx->have_accel)
    {
        matrixType ax = (matrixType) ctx->accel[0];
        matrixType ay = (matrixType) ctx->accel[1];
        matrixType az = (matrixType) ctx->accel[2];
        matrixType halfDt2 = (matrixType) 0.5 * dt * dt;
        pred[0] += halfDt2 * ax;
        pred[1] += halfDt2 * ay;
        pred[2] += halfDt2 * az;
        pred[3] += dt * ax;
        pred[4] += dt * ay;
        pred[5] += dt * az;
    }

    for (i = 0; i < 6; ++i)
    {
        SET_MATRIX(*xPredicted, i, 0, pred[i]);
    }
}

static void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    int i;
    (void) ekf;
    (void) userData;
    for (i = 0; i < 6; ++i)
    {
        SET_MATRIX(*z, i, 0, ACCESS_MATRIX(*x, i, 0));
    }
}

static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    int i;
    int j;
    (void) x;
    (void) userData;
    for (i = 0; i < 6; ++i)
    {
        for (j = 0; j < 6; ++j)
        {
            SET_MATRIX(*J, i, j, ACCESS_MATRIX(*(ekf->A), i, j));
        }
    }
}

static void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    int i;
    int j;
    (void) x;
    (void) ekf;
    (void) userData;
    for (i = 0; i < 6; ++i)
    {
        for (j = 0; j < 6; ++j)
        {
            SET_MATRIX(*J, i, j, (i == j) ? (matrixType) 1 : (matrixType) 0);
        }
    }
}

static void UpdateAMatrix(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData)
{
    int i;
    int j;
    (void) x;
    (void) ekf;
    (void) userData;
    for (i = 0; i < 6; ++i)
    {
        for (j = 0; j < 6; ++j)
        {
            SET_MATRIX(*A, i, j, (i == j) ? (matrixType) 1 : (matrixType) 0);
        }
    }
    SET_MATRIX(*A, 0, 3, (matrixType) time);
    SET_MATRIX(*A, 1, 4, (matrixType) time);
    SET_MATRIX(*A, 2, 5, (matrixType) time);
}

static void free_arrays(ImuArray* imu, GpsArray* gps, TickArray* ticks, TickArray* timeline)
{
    free(imu->values);
    free(gps->values);
    free(ticks->values);
    free(timeline->values);
    imu->values = NULL;
    gps->values = NULL;
    ticks->values = NULL;
    timeline->values = NULL;
}

static void fill_process_noise(EKFMatrix* q)
{
    int i;
    int j;
    for (i = 0; i < 6; ++i)
    {
        for (j = 0; j < 6; ++j)
        {
            SET_MATRIX(*q, i, j, (matrixType) 0);
        }
    }
    SET_MATRIX(*q, 0, 0, (matrixType) 0.05);
    SET_MATRIX(*q, 1, 1, (matrixType) 0.05);
    SET_MATRIX(*q, 2, 2, (matrixType) 0.10);
    SET_MATRIX(*q, 3, 3, (matrixType) 0.20);
    SET_MATRIX(*q, 4, 4, (matrixType) 0.20);
    SET_MATRIX(*q, 5, 5, (matrixType) 0.30);
}

static void fill_measurement_noise(EKFMatrix* r, double hdop, int fixType)
{
    int i;
    int j;
    double posScale;
    double velScale;
    double hdopClamped = hdop > 0.7 ? hdop : 0.7;
    for (i = 0; i < 6; ++i)
    {
        for (j = 0; j < 6; ++j)
        {
            SET_MATRIX(*r, i, j, (matrixType) 0);
        }
    }
    posScale = hdopClamped * hdopClamped;
    if (fixType < 3)
    {
        posScale *= 2.0;
    }
    velScale = hdopClamped * hdopClamped;
    SET_MATRIX(*r, 0, 0, (matrixType) (4.0 * posScale));
    SET_MATRIX(*r, 1, 1, (matrixType) (4.0 * posScale));
    SET_MATRIX(*r, 2, 2, (matrixType) (16.0 * posScale));
    SET_MATRIX(*r, 3, 3, (matrixType) (1.5 * velScale));
    SET_MATRIX(*r, 4, 4, (matrixType) (1.5 * velScale));
    SET_MATRIX(*r, 5, 5, (matrixType) (6.0 * velScale));
}

int main(int argc, char** argv)
{
    char imuPath[512];
    char gpsPath[512];
    const char* logsDir;
    const char* outputCsv;
    long long periodMs = DEFAULT_PERIOD_MS;
    ImuArray imu = {0};
    GpsArray gps = {0};
    TickArray allTicks = {0};
    TickArray timeline = {0};
    size_t imuIdx = 0;
    size_t gpsIdx = 0;
    bool haveImuLatest = false;
    bool haveGpsLatest = false;
    ImuRow imuLatest = {0};
    GpsRow gpsLatest = {0};
    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    EKFMeasurement measurement = {0};
    FILE* out = NULL;
    bool initialized = false;
    double refLatDeg = 0.0;
    double refLonDeg = 0.0;
    double refAltM = 0.0;
    long long prevTick = -1;
    size_t t;

    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s <logs_dir> <output_csv> [period_ms]\n", argv[0]);
        return 2;
    }

    logsDir = argv[1];
    outputCsv = argv[2];
    if (argc >= 4)
    {
        if (!parse_long_long(argv[3], &periodMs) || periodMs <= 0)
        {
            fprintf(stderr, "Invalid period_ms: %s\n", argv[3]);
            return 2;
        }
    }

    join_path(imuPath, sizeof(imuPath), logsDir, "ICM20948.CSV");
    join_path(gpsPath, sizeof(gpsPath), logsDir, "SAM_M8Q.CSV");

    if (load_imu_rows(imuPath, &imu, &allTicks) != 0)
    {
        free_arrays(&imu, &gps, &allTicks, &timeline);
        return 1;
    }
    if (load_gps_rows(gpsPath, &gps, &allTicks) != 0)
    {
        free_arrays(&imu, &gps, &allTicks, &timeline);
        return 1;
    }
    if (imu.len == 0 || gps.len == 0)
    {
        fprintf(stderr, "Required sensor rows missing. ICM=%zu GPS=%zu\n", imu.len, gps.len);
        free_arrays(&imu, &gps, &allTicks, &timeline);
        return 1;
    }

    qsort(imu.values, imu.len, sizeof(ImuRow), compare_imu_by_tick);
    qsort(gps.values, gps.len, sizeof(GpsRow), compare_gps_by_tick);
    qsort(allTicks.values, allTicks.len, sizeof(long long), compare_ll);
    if (build_timeline(&allTicks, periodMs, &timeline) != 0)
    {
        fprintf(stderr, "Failed to build timeline.\n");
        free_arrays(&imu, &gps, &allTicks, &timeline);
        return 1;
    }

    PRE_INIT_ALLOC(&ekf, 6, 6, true);
    INIT_MATRIX(options.x0, 6, 1);
    INIT_MATRIX(options.P0, 6, 6);
    INIT_MATRIX(options.Q, 6, 6);
    INIT_MATRIX(options.R, 6, 6);
    INIT_MATRIX(options.A, 6, 6);
    INIT_MATRIX(measurement.z, 6, 1);

    setIdentityMatrix(options.P0);
    fill_process_noise(options.Q);
    fill_measurement_noise(options.R, 1.0, 3);
    setIdentityMatrix(options.A);

    options.n = 6;
    options.numberOfStates = 6;
    options.numberOfMeasurements = 6;
    options.f = TransitionFunction;
    options.h = MeasurementFunction;
    options.jacobianF = StateJacobianFunction;
    options.jacobianH = MeasurementJacobianFunction;
    options.updateAMatrix = UpdateAMatrix;
    options.useFiniteDifferenceJacobian = false;
    options.mallocFlag = true;

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        fprintf(stderr, "EKFInit failed.\n");
        free_arrays(&imu, &gps, &allTicks, &timeline);
        return 1;
    }

    out = open_write_file(outputCsv);
    if (out == NULL)
    {
        fprintf(stderr, "Failed to open output CSV: %s\n", outputCsv);
        free_arrays(&imu, &gps, &allTicks, &timeline);
        return 1;
    }

    fprintf(out, "tick_ms,valid,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,ref_lat_deg,ref_lon_deg,sats,hdop\n");

    for (t = 0; t < timeline.len; ++t)
    {
        const long long tick = timeline.values[t];
        bool haveAccelNow;
        bool haveGpsNow;
        bool gpsValidNow;
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;
        double gpsPosX = 0.0;
        double gpsPosY = 0.0;
        double gpsPosZ = 0.0;
        double gpsVelX = 0.0;
        double gpsVelY = 0.0;
        double gpsVelZ = 0.0;
        int sats = 0;
        double hdop = 0.0;

        while (imuIdx < imu.len && imu.values[imuIdx].tick_ms <= tick)
        {
            imuLatest = imu.values[imuIdx];
            haveImuLatest = true;
            imuIdx++;
        }
        while (gpsIdx < gps.len && gps.values[gpsIdx].tick_ms <= tick)
        {
            gpsLatest = gps.values[gpsIdx];
            haveGpsLatest = true;
            gpsIdx++;
        }

        haveAccelNow = haveImuLatest && is_fresh(tick, imuLatest.tick_ms, IMU_FRESHNESS_MS);
        haveGpsNow = haveGpsLatest && is_fresh(tick, gpsLatest.tick_ms, GPS_FRESHNESS_MS);
        gpsValidNow = haveGpsNow && gpsLatest.valid_fix == 1;

        if (haveAccelNow)
        {
            ax = clip_value(imuLatest.accel_x_mps2, -5.0, 5.0);
            ay = clip_value(imuLatest.accel_y_mps2, -5.0, 5.0);
            az = clip_value(imuLatest.accel_z_mps2 - GRAVITY_MPS2, -5.0, 5.0);
        }

        if (haveGpsNow)
        {
            sats = gpsLatest.satellites;
            hdop = gpsLatest.hdop;
        }

        if (!initialized && gpsValidNow)
        {
            refLatDeg = gpsLatest.latitude_deg;
            refLonDeg = gpsLatest.longitude_deg;
            refAltM = gpsLatest.altitude_m;
            geodetic_to_local_m(
                gpsLatest.latitude_deg,
                gpsLatest.longitude_deg,
                gpsLatest.altitude_m,
                refLatDeg,
                refLonDeg,
                refAltM,
                &gpsPosX,
                &gpsPosY,
                &gpsPosZ);
            speed_course_to_local_velocity(gpsLatest.speed_mps, gpsLatest.course_deg, &gpsVelX, &gpsVelY, &gpsVelZ);

            SET_MATRIX(*(ekf.x), 0, 0, (matrixType) gpsPosX);
            SET_MATRIX(*(ekf.x), 1, 0, (matrixType) gpsPosY);
            SET_MATRIX(*(ekf.x), 2, 0, (matrixType) gpsPosZ);
            SET_MATRIX(*(ekf.x), 3, 0, (matrixType) gpsVelX);
            SET_MATRIX(*(ekf.x), 4, 0, (matrixType) gpsVelY);
            SET_MATRIX(*(ekf.x), 5, 0, (matrixType) gpsVelZ);

            setIdentityMatrix(ekf.P);
            {
                int i;
                for (i = 0; i < 6; ++i)
                {
                    SET_MATRIX(*(ekf.P), i, i, (matrixType) 25.0);
                }
            }
            SET_MATRIX(*(ekf.P), 0, 0, (matrixType) 100.0);
            SET_MATRIX(*(ekf.P), 1, 1, (matrixType) 100.0);
            SET_MATRIX(*(ekf.P), 2, 2, (matrixType) 100.0);
            initialized = true;
        }

        if (initialized)
        {
            double dt;
            PredictContext ctx;
            EKFReturnCodes predictRet;
            if (prevTick < 0)
            {
                dt = 0.1;
            }
            else
            {
                dt = (double) (tick - prevTick) / 1000.0;
                if (dt < 0.001)
                {
                    dt = 0.001;
                }
            }

            ctx.dt = dt;
            ctx.have_accel = haveAccelNow;
            ctx.accel[0] = ax;
            ctx.accel[1] = ay;
            ctx.accel[2] = az;
            predictRet = EKFPredict(&ekf, dt, &ctx);
            if (predictRet != EKF_SUCCESS)
            {
                fprintf(stderr, "EKFPredict failed at tick %lld\n", tick);
                fclose(out);
                free_arrays(&imu, &gps, &allTicks, &timeline);
                return 1;
            }

            if (gpsValidNow)
            {
                EKFReturnCodes updateRet;
                geodetic_to_local_m(
                    gpsLatest.latitude_deg,
                    gpsLatest.longitude_deg,
                    gpsLatest.altitude_m,
                    refLatDeg,
                    refLonDeg,
                    refAltM,
                    &gpsPosX,
                    &gpsPosY,
                    &gpsPosZ);
                speed_course_to_local_velocity(gpsLatest.speed_mps, gpsLatest.course_deg, &gpsVelX, &gpsVelY, &gpsVelZ);

                SET_MATRIX(*(measurement.z), 0, 0, (matrixType) gpsPosX);
                SET_MATRIX(*(measurement.z), 1, 0, (matrixType) gpsPosY);
                SET_MATRIX(*(measurement.z), 2, 0, (matrixType) gpsPosZ);
                SET_MATRIX(*(measurement.z), 3, 0, (matrixType) gpsVelX);
                SET_MATRIX(*(measurement.z), 4, 0, (matrixType) gpsVelY);
                SET_MATRIX(*(measurement.z), 5, 0, (matrixType) gpsVelZ);

                fill_measurement_noise(ekf.R, gpsLatest.hdop, gpsLatest.fix_type);
                updateRet = EKFUpdate(&ekf, &measurement);
                if (updateRet != EKF_SUCCESS)
                {
                    fprintf(stderr, "EKFUpdate failed at tick %lld\n", tick);
                    fclose(out);
                    free_arrays(&imu, &gps, &allTicks, &timeline);
                    return 1;
                }
            }
        }

        fprintf(
            out,
            "%lld,%d,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%d,%.10f\n",
            tick,
            initialized ? 1 : 0,
            initialized ? (double) ACCESS_MATRIX(*(ekf.x), 0, 0) : 0.0,
            initialized ? (double) ACCESS_MATRIX(*(ekf.x), 1, 0) : 0.0,
            initialized ? (double) ACCESS_MATRIX(*(ekf.x), 2, 0) : 0.0,
            initialized ? (double) ACCESS_MATRIX(*(ekf.x), 3, 0) : 0.0,
            initialized ? (double) ACCESS_MATRIX(*(ekf.x), 4, 0) : 0.0,
            initialized ? (double) ACCESS_MATRIX(*(ekf.x), 5, 0) : 0.0,
            refLatDeg,
            refLonDeg,
            sats,
            hdop);

        prevTick = tick;
    }

    fclose(out);
    EKFCleanup(&ekf);
    FREE_MATRIX(ekf.x);
    FREE_MATRIX(ekf.P);
    FREE_MATRIX(ekf.Q);
    FREE_MATRIX(ekf.R);
    FREE_MATRIX(ekf.A);
    FREE_MATRIX(options.x0);
    FREE_MATRIX(options.P0);
    FREE_MATRIX(options.Q);
    FREE_MATRIX(options.R);
    FREE_MATRIX(options.A);
    FREE_MATRIX(measurement.z);
    free_arrays(&imu, &gps, &allTicks, &timeline);

    printf("Generated EKF CSV using C EKF: %s\n", outputCsv);
    return 0;
}
