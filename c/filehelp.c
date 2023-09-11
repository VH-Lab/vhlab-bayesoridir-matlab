#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


/* lendian_fwrite - writes little endian data regardless of processor for data of 1, 2, 4 8, 16, 32, or 64 byes


*/ 

size_t lendian_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
    int x = 1;

    if ( *((char*)&x) == 1)
    {
        /* Little endian machine, use fwrite directly */
        return fwrite(ptr, size, nmemb, stream);
    }
    else
    {
        if(size == sizeof(uint8_t))     //1 Byte
        {
            return fwrite(ptr, size, nmemb, stream);
        }
        else if(size == sizeof(uint16_t))   //2 Byte
        {
            /* Big endian machine, pre-process first */
            unsigned char *buffer = malloc(size*nmemb);
            unsigned char *input = (unsigned char*) ptr;

            for (uint32_t i=0; i<nmemb; i++)
            {           
                buffer[2*i] = input[2*i + 1];
                buffer[2*i + 1] = input[2*i];
            }
            int ret =fwrite((void*)buffer, size, nmemb, stream);      
            free(buffer);
            return ret;
        }
        else if(size == sizeof(uint32_t)) //4 Byte
        { 
            /* Big endian machine, pre-process first */
            unsigned char *buffer = malloc(size*nmemb);
            unsigned char *input = (unsigned char*) ptr;

            for (uint32_t i=0; i<nmemb; i++)
            {           
                buffer[4*i    ] = input[4*i + 3];
                buffer[4*i + 1] = input[4*i + 2];
                buffer[4*i + 2] = input[4*i + 1];
                buffer[4*i + 3] = input[4*i    ];
            }

            int ret =fwrite((void*)buffer, size, nmemb, stream);      
            free(buffer);
            return ret;
        }
        else if(size == sizeof(uint64_t)) //8 Byte
        { 
            /* Big endian machine, pre-process first */
            unsigned char *buffer = malloc(size*nmemb);
            unsigned char *input = (unsigned char*) ptr;

            for (uint32_t i=0; i<nmemb; i++)
            {           
                buffer[8*i    ] = input[4*i + 7];
                buffer[8*i + 1] = input[4*i + 6];
                buffer[8*i + 2] = input[4*i + 5];
                buffer[8*i + 3] = input[4*i + 4];
                buffer[8*i + 4] = input[4*i + 3];
                buffer[8*i + 5] = input[4*i + 2];
                buffer[8*i + 6] = input[4*i + 1];
                buffer[8*i + 7] = input[4*i    ];
            }

            int ret =fwrite((void*)buffer, size, nmemb, stream);      
            free(buffer);
            return ret;
        }
        else
        {
            printf("%s Function received invalid element size:%ld\n",__FUNCTION__,size);
            return -1;
        }
 
    }  
}


size_t lendian_fread(void *ptr, size_t size, size_t nmemb,FILE *stream)
{
    int x = 1;
    int ret = 0;

    if ( *((char*)&x) == 1)
    {
        /* Little endian machine, use fwrite directly */
        return fread(ptr, size, nmemb, stream);
    }
    else
    {
        if(size == sizeof(uint8_t))     //1 Byte
        {
            return fread(ptr, size, nmemb, stream);
        }
        else if(size == sizeof(uint16_t))   //2 Byte
        {
            /* Big endian machine, pre-process first */
            unsigned char *buffer = malloc(size*nmemb);
            ret =fread((void*)buffer, size, nmemb, stream);      
            unsigned char *output = (unsigned char*) ptr;

            for (uint32_t i=0; i<nmemb; i++)
            {           
                output[2*i] = buffer[2*i + 1];
                output[2*i + 1] = buffer[2*i];
            }
            free(buffer);
            return ret;
        }
        else if(size == sizeof(uint32_t)) //4 Byte
        { 
            /* Big endian machine, pre-process first */
            unsigned char *buffer = malloc(size*nmemb);
            ret =fread((void*)buffer, size, nmemb, stream);      
            unsigned char *output = (unsigned char*) ptr;

            for (uint32_t i=0; i<nmemb; i++)
            {           
                output[4*i    ] = buffer[4*i + 3];
                output[4*i + 1] = buffer[4*i + 2];
                output[4*i + 2] = buffer[4*i + 1];
                output[4*i + 3] = buffer[4*i    ];
            }

            free(buffer);
            return ret;
        }
        else if(size == sizeof(uint64_t)) //8 Byte
        { 
            /* Big endian machine, pre-process first */

            unsigned char *buffer = malloc(size*nmemb);
            ret =fread((void*)buffer, size, nmemb, stream);      
            unsigned char *output = (unsigned char*) ptr;

            for (uint32_t i=0; i<nmemb; i++)
            {           
                output[8*i    ] = buffer[4*i + 7];
                output[8*i + 1] = buffer[4*i + 6];
                output[8*i + 2] = buffer[4*i + 5];
                output[8*i + 3] = buffer[4*i + 4];
                output[8*i + 4] = buffer[4*i + 3];
                output[8*i + 5] = buffer[4*i + 2];
                output[8*i + 6] = buffer[4*i + 1];
                output[8*i + 7] = buffer[4*i    ];
            }

            free(buffer);
            return ret;
        }
        else
        {
            printf("%s Function received invalid element size:%ld\n",__FUNCTION__,size);
            return -1;
        }
 
    }  
}



