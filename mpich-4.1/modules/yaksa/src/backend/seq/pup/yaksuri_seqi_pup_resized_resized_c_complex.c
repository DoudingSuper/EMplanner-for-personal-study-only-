/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 *
 * DO NOT EDIT: AUTOMATICALLY GENERATED FILE !!
 */

#include <string.h>
#include <stdint.h>
#include <wchar.h>
#include "yaksuri_seqi_pup.h"

int yaksuri_seqi_pack_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + i * extent)), *((float _Complex *) (void *) (dbuf + idx)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + i * extent)), *((float _Complex *) (void *) (dbuf + idx)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + i * extent)), *((float _Complex *) (void *) (dbuf + idx)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_unpack_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_pack_hvector_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    intptr_t count1 = type->u.hvector.count;
    intptr_t blocklength1 ATTRIBUTE((unused)) = type->u.hvector.blocklength;
    intptr_t stride1 = type->u.hvector.stride;
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.hvector.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.hvector.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + i * extent + j1 * stride1 + k1 * extent2)), *((float _Complex *) (void *) (dbuf + idx)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + i * extent + j1 * stride1 + k1 * extent2)), *((float _Complex *) (void *) (dbuf + idx)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + i * extent + j1 * stride1 + k1 * extent2)), *((float _Complex *) (void *) (dbuf + idx)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_unpack_hvector_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    intptr_t count1 = type->u.hvector.count;
    intptr_t blocklength1 ATTRIBUTE((unused)) = type->u.hvector.blocklength;
    intptr_t stride1 = type->u.hvector.stride;
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.hvector.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.hvector.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + j1 * stride1 + k1 * extent2)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + j1 * stride1 + k1 * extent2)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + j1 * stride1 + k1 * extent2)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_pack_blkhindx_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    intptr_t count1 = type->u.blkhindx.count;
    intptr_t blocklength1 ATTRIBUTE((unused)) = type->u.blkhindx.blocklength;
    intptr_t *restrict array_of_displs1 = type->u.blkhindx.array_of_displs;
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.blkhindx.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.blkhindx.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + i * extent + array_of_displs1[j1] + k1 * extent2)), *((float _Complex *) (void *) (dbuf + idx)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + i * extent + array_of_displs1[j1] + k1 * extent2)), *((float _Complex *) (void *) (dbuf + idx)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + i * extent + array_of_displs1[j1] + k1 * extent2)), *((float _Complex *) (void *) (dbuf + idx)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_unpack_blkhindx_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    intptr_t count1 = type->u.blkhindx.count;
    intptr_t blocklength1 ATTRIBUTE((unused)) = type->u.blkhindx.blocklength;
    intptr_t *restrict array_of_displs1 = type->u.blkhindx.array_of_displs;
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.blkhindx.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.blkhindx.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + array_of_displs1[j1] + k1 * extent2)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + array_of_displs1[j1] + k1 * extent2)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < blocklength1; k1++) {
                        YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + array_of_displs1[j1] + k1 * extent2)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_pack_hindexed_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    intptr_t count1 = type->u.hindexed.count;
    intptr_t *restrict array_of_blocklengths1 = type->u.hindexed.array_of_blocklengths;
    intptr_t *restrict array_of_displs1 = type->u.hindexed.array_of_displs;
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.hindexed.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.hindexed.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < array_of_blocklengths1[j1]; k1++) {
                        YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + i * extent + array_of_displs1[j1] + k1 * extent2)), *((float _Complex *) (void *) (dbuf + idx)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < array_of_blocklengths1[j1]; k1++) {
                        YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + i * extent + array_of_displs1[j1] + k1 * extent2)), *((float _Complex *) (void *) (dbuf + idx)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < array_of_blocklengths1[j1]; k1++) {
                        YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + i * extent + array_of_displs1[j1] + k1 * extent2)), *((float _Complex *) (void *) (dbuf + idx)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_unpack_hindexed_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    intptr_t count1 = type->u.hindexed.count;
    intptr_t *restrict array_of_blocklengths1 = type->u.hindexed.array_of_blocklengths;
    intptr_t *restrict array_of_displs1 = type->u.hindexed.array_of_displs;
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.hindexed.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.hindexed.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < array_of_blocklengths1[j1]; k1++) {
                        YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + array_of_displs1[j1] + k1 * extent2)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < array_of_blocklengths1[j1]; k1++) {
                        YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + array_of_displs1[j1] + k1 * extent2)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    for (intptr_t k1 = 0; k1 < array_of_blocklengths1[j1]; k1++) {
                        YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + array_of_displs1[j1] + k1 * extent2)));
                        idx += sizeof(float _Complex);
                    }
                }
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_pack_contig_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    intptr_t count1 = type->u.contig.count;
    intptr_t stride1 = type->u.contig.child->extent;
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.contig.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.contig.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + i * extent + j1 * stride1)), *((float _Complex *) (void *) (dbuf + idx)));
                    idx += sizeof(float _Complex);
                }
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + i * extent + j1 * stride1)), *((float _Complex *) (void *) (dbuf + idx)));
                    idx += sizeof(float _Complex);
                }
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + i * extent + j1 * stride1)), *((float _Complex *) (void *) (dbuf + idx)));
                    idx += sizeof(float _Complex);
                }
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_unpack_contig_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    intptr_t count1 = type->u.contig.count;
    intptr_t stride1 = type->u.contig.child->extent;
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.contig.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.contig.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + j1 * stride1)));
                    idx += sizeof(float _Complex);
                }
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + j1 * stride1)));
                    idx += sizeof(float _Complex);
                }
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                for (intptr_t j1 = 0; j1 < count1; j1++) {
                    YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent + j1 * stride1)));
                    idx += sizeof(float _Complex);
                }
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_pack_resized_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.resized.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.resized.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + i * extent)), *((float _Complex *) (void *) (dbuf + idx)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + i * extent)), *((float _Complex *) (void *) (dbuf + idx)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + i * extent)), *((float _Complex *) (void *) (dbuf + idx)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

int yaksuri_seqi_unpack_resized_resized_resized_c_complex(const void *inbuf, void *outbuf, uintptr_t count, yaksi_type_s * type, yaksa_op_t op)
{
    int rc = YAKSA_SUCCESS;
    const char *restrict sbuf = (const char *) inbuf;
    char *restrict dbuf = (char *) outbuf;
    uintptr_t extent ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent1 ATTRIBUTE((unused)) = type->extent;
    
    uintptr_t extent2 ATTRIBUTE((unused)) = type->u.resized.child->extent;
    
    uintptr_t extent3 ATTRIBUTE((unused)) = type->u.resized.child->u.resized.child->extent;
    
    uintptr_t idx = 0;
    switch (op) {
        case YAKSA_OP__SUM:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_SUM(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        case YAKSA_OP__PROD:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_PROD(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        case YAKSA_OP__REPLACE:
        {
            for (intptr_t i = 0; i < count; i++) {
                YAKSURI_SEQI_OP_REPLACE(*((const float _Complex *) (const void *) (sbuf + idx)), *((float _Complex *) (void *) (dbuf + i * extent)));
                idx += sizeof(float _Complex);
            }
            break;
        }
        
        default:
            break;
    }
    
    return rc;
}

