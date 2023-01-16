
#if WIN32 || WIN64

    static inline FILE* FOPEN(char const* _FileName, char const* _Mode)
    {
        FILE* ret = NULL;
        fopen_s(&ret, _FileName, _Mode);
        return ret;
    }

#define STRCPY strcpy_s
#define STRCAT strcat_s
#define SPRINTF sprintf_s
#define FSCANF fscanf_s

#else
#define FOPEN fopen
#define STRCPY(dest, length, source) strcpy(dest, source)
#define STRCAT(dest, length, source) strcat(dest, source)
#define SPRINTF(dest, length, format, a) sprintf(dest, format, a)
#define SPRINTF(dest, length, format, a, b) sprintf(dest, format, a, b)
#define SPRINTF(dest, length, format, a, b, c) sprintf(dest, format, a, b, c)
#define SPRINTF(dest, length, format, a, b, c, d) sprintf(dest, format, a, b, c, d)
#define SPRINTF(dest, length, format, a, b, c, d, e) sprintf(dest, format, a, b, c, d, e)
#define FSCANF fscanf
#endif
