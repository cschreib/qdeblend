IMG=${1}
PSF=${2}
OVERSAMPLE=${3}
ERROR=${4}
XSHIFT=${5}
YSHIFT=${6}
RAC=${7}
DECC=${8}
LARGE_SCALE_FILTER=${9}
MASK_MORE=${10}
GET_SHIFT=${11}
SHIFT_ANCHOR=${12}
SHIFT_MASK=${13}
SHIFT_STEP=${14}
SHIFT_MAX=${15}

echo ${IMG}

DIR=`dirname ${IMG}`

if [ ${LARGE_SCALE_FILTER} -gt 0 ]; then
    # Backup the original image so we don't apply twice
    if [ -f ${IMG}-orig.fits ]; then
        cp ${IMG}-orig.fits ${IMG}.fits
    else
        cp ${IMG}.fits ${IMG}-orig.fits
    fi

    ovif ./make_mask ${IMG}.fits ${DIR}/mask_all.reg ${IMG}-mask-all.fits
    if [ -f ${IMG}-bad.reg ]; then
        ovif ./add_mask ${IMG}-mask-all.fits ${IMG}-bad.reg
    fi

    ovif ./large_scale_filter image=${IMG}.fits size=${LARGE_SCALE_FILTER} mask=${IMG}-mask-all.fits
fi

if [ ${MASK_MORE} -eq 0 ]; then
    ovif ./make_mask ${IMG}.fits ${DIR}/mask_border.reg ${IMG}-mask-border.fits
fi
if [ ${MASK_MORE} -eq 1 ]; then
    ovif ./make_mask ${IMG}.fits ${DIR}/mask_border_more.reg ${IMG}-mask-border.fits
fi

if [ -f ${IMG}-bad.reg ]; then
    ovif ./add_mask ${IMG}-mask-border.fits ${IMG}-bad.reg
fi

if [ ${GET_SHIFT} -eq 0 ]; then
    ovif ./multi_subtract image=${IMG}.fits psf=${DIR}/${PSF} oversample=${OVERSAMPLE} \
        mask=${IMG}-mask-border.fits out=${IMG}-allsub.fits error=${ERROR} \
        model_dir=${DIR}/imfit/models xshift=${XSHIFT} yshift=${YSHIFT} reuse_models=0 \
        ra_center=${RAC} dec_center=${DECC}
else
    ovif ./multi_subtract image=${IMG}.fits psf=${DIR}/${PSF} oversample=${OVERSAMPLE} \
        mask=${IMG}-mask-border.fits out=${IMG}-allsub.fits error=${ERROR} \
        model_dir=${DIR}/imfit/models shift_anchors=${SHIFT_ANCHOR} \
        shift_fit_radius=${SHIFT_MASK} max_shift=${SHIFT_MAX} shift_dp=${SHIFT_STEP} \
        ra_center=${RAC} dec_center=${DECC}
fi
