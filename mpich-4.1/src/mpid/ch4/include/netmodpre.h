/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef NETMODPRE_H_INCLUDED
#define NETMODPRE_H_INCLUDED

/* *INDENT-OFF* */
#include "../netmod/ofi/ofi_pre.h"
/* *INDENT-ON* */

#define MPIDI_NM_REQUEST_AM_DECL MPIDI_OFI_am_request_t ofi;
#define MPIDI_NM_REQUEST_DECL    MPIDI_OFI_request_t ofi;

#define MPIDI_NM_COMM_DECL       MPIDI_OFI_comm_t ofi;
#define MPIDI_NM_DT_DECL         MPIDI_OFI_dt_t ofi;
#define MPIDI_NM_WIN_DECL        MPIDI_OFI_win_t ofi;
#define MPIDI_NM_ADDR_DECL    MPIDI_OFI_addr_t ofi;
#define MPIDI_NM_OP_DECL         MPIDI_OFI_op_t ofi;
#define MPIDI_NM_PART_DECL         MPIDI_OFI_part_t ofi;

#endif /* NETMODPRE_H_INCLUDED */
