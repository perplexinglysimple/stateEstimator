#include "EKF.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#define BENCH_STATE_DIM 6
#define BENCH_MEAS_DIM 3
#define CALIBRATION_SAMPLES 32
#define LINE_BUFFER_SIZE 512

typedef struct ImuRow_
{
    long tick_ms;
    double accel_x_mps2;
    double accel_y_mps2;
    double accel_z_mps2;
} ImuRow;

typedef struct BmpRow_
{
    long tick_ms;
    double pressure_pa;
} BmpRow;

typedef struct GpsRow_
{
    long tick_ms;
    int valid_fix;
    double latitude_deg;
    double longitude_deg;
    double altitude_m;
} GpsRow;

typedef struct MotionInput_
{
    double ax;
    double ay;
    double az;
    double dt;
} MotionInput;

static int readNextImu(FILE* file, ImuRow* row);
static int readNextBmp(FILE* file, BmpRow* row);
static int readNextGps(FILE* file, GpsRow* row);
static void geodeticToLocalMeters(double lat_deg, double lon_deg, double lat0_deg, double lon0_deg, double* east_m,
                                  double* north_m);
static double pressureToRelativeAltitude(double pressure_pa, double pressure0_pa);
static void TransitionFunction(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData);
static void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData);
static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);
static void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData);
static void UpdateAMatrix(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData);

int main(int argc, char** argv)
{
    const char* outPath = "ekf_bench_compare.csv";
    const char* imuPath = "test_data/ICM20948.CSV";
    const char* gpsPath = "test_data/SAM_M8Q.CSV";
    const char* bmpPath = "test_data/BMP384.CSV";
    if (argc > 1 && argv[1] != NULL)
    {
        outPath = argv[1];
    }
    if (argc > 2 && argv[2] != NULL)
    {
        imuPath = argv[2];
    }
    if (argc > 3 && argv[3] != NULL)
    {
        gpsPath = argv[3];
    }
    if (argc > 4 && argv[4] != NULL)
    {
        bmpPath = argv[4];
    }

    FILE* imuFile = fopen(imuPath, "r");
    FILE* gpsFile = fopen(gpsPath, "r");
    FILE* bmpFile = fopen(bmpPath, "r");
    FILE* outFile = fopen(outPath, "w");
    if (!imuFile || !gpsFile || !bmpFile || !outFile)
    {
        LOG_ERROR("Failed to open one of the bench-data files.");
        if (imuFile)
        {
            fclose(imuFile);
        }
        if (gpsFile)
        {
            fclose(gpsFile);
        }
        if (bmpFile)
        {
            fclose(bmpFile);
        }
        if (outFile)
        {
            fclose(outFile);
        }
        return -1;
    }

    char header[LINE_BUFFER_SIZE];
    if (!fgets(header, sizeof(header), imuFile) || !fgets(header, sizeof(header), gpsFile) ||
        !fgets(header, sizeof(header), bmpFile))
    {
        LOG_ERROR("Failed to read CSV headers.");
        fclose(imuFile);
        fclose(gpsFile);
        fclose(bmpFile);
        fclose(outFile);
        return -1;
    }

    EKFState ekf = {0};
    EKFConfigOptions options = {0};
    PRE_INIT_ALLOC(&ekf, BENCH_STATE_DIM, BENCH_MEAS_DIM, true);

    STATIC_MATRIX_DIRECTIVE(options.x0, BENCH_STATE_DIM, 1, x0);
    STATIC_MATRIX_DIRECTIVE(options.P0, BENCH_STATE_DIM, BENCH_STATE_DIM, P0);
    STATIC_MATRIX_DIRECTIVE(options.Q, BENCH_STATE_DIM, BENCH_STATE_DIM, Q);
    STATIC_MATRIX_DIRECTIVE(options.R, BENCH_MEAS_DIM, BENCH_MEAS_DIM, R);
    STATIC_MATRIX_DIRECTIVE(options.A, BENCH_STATE_DIM, BENCH_STATE_DIM, A);

    for (int i = 0; i < BENCH_STATE_DIM; ++i)
    {
        SET_MATRIX(*(options.P0), i, i, (i < 3) ? 25.0 : 4.0);
        SET_MATRIX(*(options.Q), i, i, (i < 3) ? 0.05 : 0.2);
        SET_MATRIX(*(options.A), i, i, 1.0);
    }
    SET_MATRIX(*(options.R), 0, 0, 9.0);
    SET_MATRIX(*(options.R), 1, 1, 9.0);
    SET_MATRIX(*(options.R), 2, 2, 4.0);

    options.n = BENCH_STATE_DIM;
    options.f = TransitionFunction;
    options.h = MeasurementFunction;
    options.updateAMatrix = UpdateAMatrix;
    options.numberOfStates = BENCH_STATE_DIM;
    options.numberOfMeasurements = BENCH_MEAS_DIM;
    options.useFiniteDifferenceJacobian = false;
    options.jacobianF = StateJacobianFunction;
    options.jacobianH = MeasurementJacobianFunction;
    options.mallocFlag = true;

    if (EKFInit(&ekf, &options) != EKF_SUCCESS)
    {
        LOG_ERROR("EKFInit() failed in bench export test.");
        fclose(imuFile);
        fclose(gpsFile);
        fclose(bmpFile);
        fclose(outFile);
        return -1;
    }

    EKFMeasurement measurement = {0};
    STATIC_MATRIX_DIRECTIVE(measurement.z, BENCH_MEAS_DIM, 1, z);

    fprintf(outFile, "step,tick_ms,x0,x1,x2,x3,x4,x5,P0,P1,P2,P3,P4,P5\n");

    ImuRow imuRow = {0};
    BmpRow bmpRow = {0};
    GpsRow gpsRow = {0};
    int haveBmp = readNextBmp(bmpFile, &bmpRow);
    int haveGps = readNextGps(gpsFile, &gpsRow);

    double accelBias[3] = {0.0, 0.0, 0.0};
    int biasSamples = 0;
    long previousTick = -1;

    double pressure0 = 0.0;
    int havePressure0 = 0;
    double latestBaroZ = 0.0;
    int haveBaro = 0;

    double originLatDeg = 0.0;
    double originLonDeg = 0.0;
    int haveGpsOrigin = 0;
    double latestGpsX = 0.0;
    double latestGpsY = 0.0;
    int haveGpsMeasurement = 0;

    int measurementDirty = 0;
    int step = 0;
    while (readNextImu(imuFile, &imuRow))
    {
        while (haveBmp && bmpRow.tick_ms <= imuRow.tick_ms)
        {
            if (!havePressure0)
            {
                pressure0 = bmpRow.pressure_pa;
                havePressure0 = 1;
            }
            latestBaroZ = pressureToRelativeAltitude(bmpRow.pressure_pa, pressure0);
            haveBaro = 1;
            measurementDirty = 1;
            haveBmp = readNextBmp(bmpFile, &bmpRow);
        }

        while (haveGps && gpsRow.tick_ms <= imuRow.tick_ms)
        {
            if (gpsRow.valid_fix && gpsRow.latitude_deg != 0.0 && gpsRow.longitude_deg != 0.0)
            {
                if (!haveGpsOrigin)
                {
                    originLatDeg = gpsRow.latitude_deg;
                    originLonDeg = gpsRow.longitude_deg;
                    haveGpsOrigin = 1;
                }
                geodeticToLocalMeters(gpsRow.latitude_deg, gpsRow.longitude_deg, originLatDeg, originLonDeg,
                                      &latestGpsX, &latestGpsY);
                haveGpsMeasurement = 1;
                measurementDirty = 1;
            }
            haveGps = readNextGps(gpsFile, &gpsRow);
        }

        if (biasSamples < CALIBRATION_SAMPLES)
        {
            accelBias[0] += imuRow.accel_x_mps2;
            accelBias[1] += imuRow.accel_y_mps2;
            accelBias[2] += imuRow.accel_z_mps2;
            ++biasSamples;
        }

        if (previousTick < 0)
        {
            previousTick = imuRow.tick_ms;
            continue;
        }

        double dt = (double) (imuRow.tick_ms - previousTick) / 1000.0;
        previousTick = imuRow.tick_ms;
        if (dt <= 0.0)
        {
            continue;
        }

        double biasScale = (biasSamples > 0) ? 1.0 / (double) biasSamples : 0.0;
        MotionInput input = {
            .ax = imuRow.accel_x_mps2 - accelBias[0] * biasScale,
            .ay = imuRow.accel_y_mps2 - accelBias[1] * biasScale,
            .az = imuRow.accel_z_mps2 - accelBias[2] * biasScale,
            .dt = dt,
        };

        if (EKFPredict(&ekf, dt, &input) != EKF_SUCCESS)
        {
            LOG_ERROR("EKFPredict() failed in bench export test.");
            EKFCleanup(&ekf);
            fclose(imuFile);
            fclose(gpsFile);
            fclose(bmpFile);
            fclose(outFile);
            return -1;
        }

        if (measurementDirty && haveBaro && haveGpsMeasurement)
        {
            SET_MATRIX(*(measurement.z), 0, 0, latestGpsX);
            SET_MATRIX(*(measurement.z), 1, 0, latestGpsY);
            SET_MATRIX(*(measurement.z), 2, 0, latestBaroZ);
            if (EKFUpdate(&ekf, &measurement) != EKF_SUCCESS)
            {
                LOG_ERROR("EKFUpdate() failed in bench export test.");
                EKFCleanup(&ekf);
                fclose(imuFile);
                fclose(gpsFile);
                fclose(bmpFile);
                fclose(outFile);
                return -1;
            }
            measurementDirty = 0;
        }

        fprintf(outFile, "%d,%ld", step, imuRow.tick_ms);
        for (int i = 0; i < BENCH_STATE_DIM; ++i)
        {
            fprintf(outFile, ",%.10f", ACCESS_MATRIX(*(ekf.x), i, 0));
        }
        for (int i = 0; i < BENCH_STATE_DIM; ++i)
        {
            fprintf(outFile, ",%.10f", ACCESS_MATRIX(*(ekf.P), i, i));
        }
        fprintf(outFile, "\n");
        ++step;
    }

    EKFCleanup(&ekf);
    fclose(imuFile);
    fclose(gpsFile);
    fclose(bmpFile);
    fclose(outFile);
    return 0;
}

static int readNextImu(FILE* file, ImuRow* row)
{
    char line[LINE_BUFFER_SIZE];
    while (fgets(line, sizeof(line), file))
    {
        long tick_ms = 0;
        double ax = 0.0;
        double ay = 0.0;
        double az = 0.0;
        if (sscanf(line, "%ld,%lf,%lf,%lf", &tick_ms, &ax, &ay, &az) == 4)
        {
            row->tick_ms = tick_ms;
            row->accel_x_mps2 = ax;
            row->accel_y_mps2 = ay;
            row->accel_z_mps2 = az;
            return 1;
        }
    }
    return 0;
}

static int readNextBmp(FILE* file, BmpRow* row)
{
    char line[LINE_BUFFER_SIZE];
    while (fgets(line, sizeof(line), file))
    {
        long tick_ms = 0;
        double temperature_c = 0.0;
        double pressure_pa = 0.0;
        if (sscanf(line, "%ld,%lf,%lf", &tick_ms, &temperature_c, &pressure_pa) == 3)
        {
            row->tick_ms = tick_ms;
            row->pressure_pa = pressure_pa;
            return 1;
        }
    }
    return 0;
}

static int readNextGps(FILE* file, GpsRow* row)
{
    char line[LINE_BUFFER_SIZE];
    while (fgets(line, sizeof(line), file))
    {
        long tick_ms = 0;
        int valid_fix = 0;
        int fix_type = 0;
        int satellites = 0;
        double latitude_deg = 0.0;
        double longitude_deg = 0.0;
        double altitude_m = 0.0;
        double speed_mps = 0.0;
        double course_deg = 0.0;
        double hdop = 0.0;
        if (sscanf(line, "%ld,%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf", &tick_ms, &valid_fix, &fix_type, &satellites,
                   &latitude_deg, &longitude_deg, &altitude_m, &speed_mps, &course_deg, &hdop) == 10)
        {
            (void) fix_type;
            (void) satellites;
            (void) altitude_m;
            (void) speed_mps;
            (void) course_deg;
            (void) hdop;
            row->tick_ms = tick_ms;
            row->valid_fix = valid_fix;
            row->latitude_deg = latitude_deg;
            row->longitude_deg = longitude_deg;
            row->altitude_m = altitude_m;
            return 1;
        }
    }
    return 0;
}

static void geodeticToLocalMeters(double lat_deg, double lon_deg, double lat0_deg, double lon0_deg, double* east_m,
                                  double* north_m)
{
    const double earthRadiusM = 6378137.0;
    const double degToRad = 3.14159265358979323846 / 180.0;
    double lat = lat_deg * degToRad;
    double lon = lon_deg * degToRad;
    double lat0 = lat0_deg * degToRad;
    double lon0 = lon0_deg * degToRad;
    *east_m = earthRadiusM * cos(lat0) * (lon - lon0);
    *north_m = earthRadiusM * (lat - lat0);
}

static double pressureToRelativeAltitude(double pressure_pa, double pressure0_pa)
{
    if (pressure_pa <= 0.0 || pressure0_pa <= 0.0)
    {
        return 0.0;
    }
    return 44330.0 * (1.0 - pow(pressure_pa / pressure0_pa, 0.19029495718363465));
}

static void TransitionFunction(EKFMatrix* x, EKFMatrix* x_pred, EKFState* ekf, void* userData)
{
    (void) ekf;
    const MotionInput* input = (const MotionInput*) userData;
    const double dt = input ? input->dt : 0.0;
    const double ax = input ? input->ax : 0.0;
    const double ay = input ? input->ay : 0.0;
    const double az = input ? input->az : 0.0;

    double px = ACCESS_MATRIX(*x, 0, 0);
    double py = ACCESS_MATRIX(*x, 1, 0);
    double pz = ACCESS_MATRIX(*x, 2, 0);
    double vx = ACCESS_MATRIX(*x, 3, 0);
    double vy = ACCESS_MATRIX(*x, 4, 0);
    double vz = ACCESS_MATRIX(*x, 5, 0);

    SET_MATRIX(*x_pred, 0, 0, px + vx * dt + 0.5 * ax * dt * dt);
    SET_MATRIX(*x_pred, 1, 0, py + vy * dt + 0.5 * ay * dt * dt);
    SET_MATRIX(*x_pred, 2, 0, pz + vz * dt + 0.5 * az * dt * dt);
    SET_MATRIX(*x_pred, 3, 0, vx + ax * dt);
    SET_MATRIX(*x_pred, 4, 0, vy + ay * dt);
    SET_MATRIX(*x_pred, 5, 0, vz + az * dt);
}

static void MeasurementFunction(EKFMatrix* x, EKFMatrix* z, EKFState* ekf, void* userData)
{
    (void) ekf;
    (void) userData;
    for (int i = 0; i < BENCH_MEAS_DIM; ++i)
    {
        SET_MATRIX(*z, i, 0, ACCESS_MATRIX(*x, i, 0));
    }
}

static void StateJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void) x;
    (void) userData;
    for (int i = 0; i < BENCH_STATE_DIM; ++i)
    {
        for (int j = 0; j < BENCH_STATE_DIM; ++j)
        {
            SET_MATRIX(*J, i, j, ACCESS_MATRIX(*(ekf->A), i, j));
        }
    }
}

static void MeasurementJacobianFunction(EKFMatrix* x, EKFMatrix* J, EKFState* ekf, void* userData)
{
    (void) x;
    (void) ekf;
    (void) userData;
    for (int i = 0; i < BENCH_MEAS_DIM; ++i)
    {
        for (int j = 0; j < BENCH_STATE_DIM; ++j)
        {
            SET_MATRIX(*J, i, j, (i == j) ? 1.0 : 0.0);
        }
    }
}

static void UpdateAMatrix(EKFMatrix* A, EKFMatrix* x, struct EKFState_* ekf, double time, void* userData)
{
    (void) x;
    (void) ekf;
    (void) userData;
    for (int i = 0; i < BENCH_STATE_DIM; ++i)
    {
        for (int j = 0; j < BENCH_STATE_DIM; ++j)
        {
            SET_MATRIX(*A, i, j, 0.0);
        }
    }
    for (int i = 0; i < BENCH_STATE_DIM; ++i)
    {
        SET_MATRIX(*A, i, i, 1.0);
    }
    SET_MATRIX(*A, 0, 3, time);
    SET_MATRIX(*A, 1, 4, time);
    SET_MATRIX(*A, 2, 5, time);
}
