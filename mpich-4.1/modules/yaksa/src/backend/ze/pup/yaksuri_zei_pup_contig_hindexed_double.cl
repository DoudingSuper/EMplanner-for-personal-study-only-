/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 *
 * DO NOT EDIT: AUTOMATICALLY GENERATED FILE !!
 */

typedef signed char int8_t;
typedef signed short int int16_t;
typedef signed int int32_t;
typedef signed long int64_t;
typedef unsigned char uint8_t;
typedef unsigned short int uint16_t;
typedef unsigned int uint32_t;
typedef unsigned long uint64_t;
#include "yaksuri_zei_md.h"

__kernel void yaksuri_zei_kernel_pack_SUM_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) += *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_SUM_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double))) += *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_PROD_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) *= *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_PROD_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double))) *= *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_REPLACE_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) = *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_REPLACE_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double))) = *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_MIN_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double))) > *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MIN_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) > *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_MAX_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double))) < *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MAX_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t *array_of_displs2 = md->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) < *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + array_of_displs2[x2] + x3 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_SUM_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) += *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_SUM_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) += *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_PROD_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) *= *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_PROD_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) *= *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_REPLACE_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) = *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_REPLACE_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) = *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_MIN_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) > *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MIN_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) > *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_MAX_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) < *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MAX_hvector_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hvector.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hvector.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hvector.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hvector.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t stride1 = md->u.hvector.stride;
    intptr_t stride2 = md->u.hvector.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hvector.child->extent;
    intptr_t *array_of_displs3 = md->u.hvector.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) < *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_SUM_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) += *((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_SUM_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) += *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_PROD_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) *= *((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_PROD_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) *= *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_REPLACE_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) = *((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_REPLACE_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) = *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_MIN_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) > *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MIN_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) > *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_MAX_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) < *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MAX_blkhindx_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.blocklength;
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.blkhindx.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.blkhindx.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.blkhindx.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.blkhindx.array_of_displs;
    intptr_t stride2 = md->u.blkhindx.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.blkhindx.child->extent;
    intptr_t *array_of_displs3 = md->u.blkhindx.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) < *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_SUM_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) += *((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_SUM_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) += *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_PROD_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) *= *((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_PROD_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) *= *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_REPLACE_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) = *((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_REPLACE_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) = *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_MIN_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) > *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MIN_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) > *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_MAX_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) < *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MAX_hindexed_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.count;
    
    uintptr_t x1;
    for (intptr_t i = 0; i < md->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x1 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.count;
    
    uintptr_t x3 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.hindexed.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x4;
    for (intptr_t i = 0; i < md->u.hindexed.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x4 = i;
                    res %= in_elems;
                    inner_elements = md->u.hindexed.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x5 = res;
    
    intptr_t *array_of_displs1 = md->u.hindexed.array_of_displs;
    intptr_t stride2 = md->u.hindexed.child->u.contig.child->extent;
    uintptr_t extent2 = md->u.hindexed.child->extent;
    intptr_t *array_of_displs3 = md->u.hindexed.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) < *((double *) (void *) (dbuf + x0 * extent + array_of_displs1[x1] + x2 * extent2 + x3 * stride2 + array_of_displs3[x4] + x5 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_SUM_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) += *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_SUM_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double))) += *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_PROD_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) *= *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_PROD_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double))) *= *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_REPLACE_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) = *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_REPLACE_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double))) = *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_MIN_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double))) > *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MIN_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) > *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_MAX_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double))) < *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MAX_contig_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.count;
    
    uintptr_t x2 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.contig.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x3;
    for (intptr_t i = 0; i < md->u.contig.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.contig.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x3 = i;
                    res %= in_elems;
                    inner_elements = md->u.contig.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x4 = res;
    
    intptr_t stride1 = md->u.contig.child->extent;
    intptr_t stride2 = md->u.contig.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.contig.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) < *((double *) (void *) (dbuf + x0 * extent + x1 * stride1 + x2 * stride2 + array_of_displs3[x3] + x4 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_SUM_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) += *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_SUM_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double))) += *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_PROD_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) *= *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_PROD_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double))) *= *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_REPLACE_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + idx * sizeof(double))) = *((const double *) (const void *) (sbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_unpack_REPLACE_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double))) = *((const double *) (const void *) (sbuf + idx * sizeof(double)));
}

__kernel void yaksuri_zei_kernel_pack_MIN_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double))) > *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MIN_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) > *((double *) (void *) (dbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_pack_MAX_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    sbuf = (__global const char *) ((__global char *)sbuf - md->true_lb);
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double))), *((double *) (void *) (dbuf + idx * sizeof(double)))};
    *((double *) (void *) (dbuf + idx * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double))) < *((double *) (void *) (dbuf + idx * sizeof(double)))];
}

__kernel void yaksuri_zei_kernel_unpack_MAX_resized_contig_hindexed_double(__global const void *inbuf, __global void *outbuf, unsigned long count, __global const yaksuri_zei_md_s *__restrict__ md)
{
    __global const char *__restrict__ sbuf = (__global char *) inbuf;
    __global char *__restrict__ dbuf = (__global char *) outbuf;
    dbuf = dbuf - md->true_lb;
    uintptr_t extent = md->extent;
    uintptr_t idx = get_global_id(0);
    uintptr_t res = idx;
    uintptr_t inner_elements = md->num_elements;
    
    if (idx >= (count * inner_elements))
        return;
    
    uintptr_t x0 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.count;
    
    uintptr_t x1 = res / inner_elements;
    res %= inner_elements;
    inner_elements /= md->u.resized.child->u.contig.child->u.hindexed.count;
    
    uintptr_t x2;
    for (intptr_t i = 0; i < md->u.resized.child->u.contig.child->u.hindexed.count; i++) {
            uintptr_t in_elems = md->u.resized.child->u.contig.child->u.hindexed.array_of_blocklengths[i] *
                                 md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
            if (res < in_elems) {
                    x2 = i;
                    res %= in_elems;
                    inner_elements = md->u.resized.child->u.contig.child->u.hindexed.child->num_elements;
                    break;
            } else {
                    res -= in_elems;
            }
    }
    
    uintptr_t x3 = res;
    
    intptr_t stride2 = md->u.resized.child->u.contig.child->extent;
    intptr_t *array_of_displs3 = md->u.resized.child->u.contig.child->u.hindexed.array_of_displs;
    double x_[2] = {*((const double *) (const void *) (sbuf + idx * sizeof(double))), *((double *) (void *) (dbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double)))};
    *((double *) (void *) (dbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double))) = x_[*((const double *) (const void *) (sbuf + idx * sizeof(double))) < *((double *) (void *) (dbuf + x0 * extent + x1 * stride2 + array_of_displs3[x2] + x3 * sizeof(double)))];
}

