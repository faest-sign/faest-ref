#ifndef VK_BOX_H
#define VK_BOX_H

#include "fields.h"
#include "parameters.h"

bf128_t* get_vk_128(vbb_t* vbb, unsigned int idx);
bf192_t* get_vk_192(vbb_t* vbb, unsigned int idx);
bf256_t* get_vk_256(vbb_t* vbb, unsigned int idx);

#endif