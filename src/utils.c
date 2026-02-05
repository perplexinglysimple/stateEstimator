#include "utils.h"

#ifdef __linux__
const char* GetFileNameFromPath(const char* path)
{
    // basename may modify its input, so copy to a static buffer.
    static char pathCopy[4096];
    if (path == NULL)
    {
        return NULL;
    }
    strncpy(pathCopy, path, sizeof(pathCopy) - 1);
    pathCopy[sizeof(pathCopy) - 1] = '\0';
    return basename(pathCopy);
}
#elif _WIN32

// This is stupid windows stuff to get the filename from the path.... Linux is so much simpler

#define FAILURE_NULL_ARGUMENT       ((DWORD)-1)
#define FAILURE_API_CALL            ((DWORD)-2)
#define FAILURE_INSUFFICIENT_BUFFER ((DWORD)-3)

DWORD GetBasePathFromPathName(LPCTSTR szPathName, LPTSTR szBasePath, DWORD dwBasePathSize)
{
  TCHAR   szDrive[_MAX_DRIVE] = { 0 };
  TCHAR   szDir[_MAX_DIR]     = { 0 };
  TCHAR   szFname[_MAX_FNAME] = { 0 };
  TCHAR   szExt[_MAX_EXT]     = { 0 };
  size_t  PathLength;
  DWORD   dwReturnCode;
  DWORD   dwFileNameIndex;

  // Parameter validation
  if (szPathName == NULL || szBasePath == NULL)
  {
    return FAILURE_NULL_ARGUMENT;
  }

  // Split the path into its components
  dwReturnCode = _tsplitpath_s(szPathName, szDrive, _MAX_DRIVE, szDir, _MAX_DIR, szFname, _MAX_FNAME, szExt, _MAX_EXT);
  if (dwReturnCode != 0)
  {
    _ftprintf(stderr, TEXT("Error splitting path. _tsplitpath_s returned %d.\n"), dwReturnCode);
    return FAILURE_API_CALL;
  }

  // Check that the provided buffer is large enough to store the results and a terminal null character
  PathLength = _tcslen(szDrive) + _tcslen(szDir);
  if ((PathLength + sizeof(TCHAR)) > dwBasePathSize)
  {
    _ftprintf(stderr, TEXT("Insufficient buffer. Required %zd. Provided: %d\n"), PathLength, dwBasePathSize);
    return FAILURE_INSUFFICIENT_BUFFER;
  }

  // Copy the szDrive and szDir into the provided buffer to form the basepath
  if ((dwReturnCode = _tcscpy_s(szBasePath, dwBasePathSize, szDrive)) != 0)
  {
    _ftprintf(stderr, TEXT("Error copying string. _tcscpy_s returned %d\n"), dwReturnCode);
    return FAILURE_API_CALL;
  }
  if ((dwReturnCode = _tcscat_s(szBasePath, dwBasePathSize, szDir)) != 0)
  {
    _ftprintf(stderr, TEXT("Error copying string. _tcscat_s returned %d\n"), dwReturnCode);
    return FAILURE_API_CALL;
  }

  // Calculate the index of the start of the filename
  dwFileNameIndex = (DWORD) _tcslen(szBasePath);

  return dwFileNameIndex;
}

const char* GetFileNameFromPath(const char* path)
{
    char buffer[MAX_PATH_LENGTH];
    int filenamePos = 0;

    if ((filenamePos = GetBasePathFromPathName(path, buffer, MAX_PATH_LENGTH)) < 0)
    {
        return path; // Fallback to original path if retrieval fails
    }

    return path + filenamePos;
}
#endif
