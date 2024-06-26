#include "psmx3.h"
#include "psm3_src_chksum.h"

#ifndef PSMX3_IEFS_VERSION
#define PSMX3_IEFS_VERSION	"3_2_0_0"
#endif

#ifndef PSMX3_BUILD_TIMESTAMP
#define PSMX3_BUILD_TIMESTAMP	"<Unknown>"
#endif

#ifndef PSMX3_SRC_CHECKSUM
#define PSMX3_SRC_CHECKSUM	"<Unknown>"
#endif

#ifndef PSMX3_GIT_CHECKSUM
#define PSMX3_GIT_CHECKSUM	"<Unknown>"
#endif

char psm3_IEFS_version[] = PSMX3_IEFS_VERSION;
char psm3_build_timestamp[] = PSMX3_BUILD_TIMESTAMP;
char psm3_sources_checksum[] = PSMX3_SRC_CHECKSUM;
char psm3_git_checksum[] = PSMX3_GIT_CHECKSUM;

#define PSM3_PROV_VER_MAJOR 3
#define PSM3_PROV_VER_MINOR 2
#define PSM3_PROV_VER_MAINT 0
#define PSM3_PROV_VER_PATCH 0

/* Leave last digit open for special use */
#define PSM3_PROV_VER(major, minor, maint, patch) \
	( ( ( ((major) * 100) + (minor)) << 16)	| ( ( ((maint) * 1000) + ((patch) * 10)) & 0xFFFF) )


static uint32_t psm3_provider_version =
	PSM3_PROV_VER(PSM3_PROV_VER_MAJOR, PSM3_PROV_VER_MINOR, PSM3_PROV_VER_MAINT, PSM3_PROV_VER_PATCH);
uint32_t get_psm3_provider_version() {
	return psm3_provider_version;
}
