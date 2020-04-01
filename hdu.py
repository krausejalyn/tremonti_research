from astropy.io import fits

hdu = fits.open('zooinverse_summary_v0.fit')
plateifu = hdu[1].data['PLATEIFU']
mangaid = hdu[1].data['MANGAID']
ra = hdu[1].data['RA']
dec = hdu[1].data['DEC']
z = hdu[1].data['Z']
fnugriz_absmag = hdu[1].data['FNUGRIZ_ABSMAG'] # array of 7 colors in this order:  FUV, NUV, u, g, r, i, z
fuv_r = fnugriz_absmag[:,0] - fnugriz_absmag[:,4] # FUV-r color 
log_mass = hdu[1].data['LOG_MASS']
subject_id = hdu[1].data['SUBJECT_ID']
nclass = hdu[1].data['NCLASS']
bad_re = hdu[1].data['BAD_RE']
bad_re_err = hdu[1].data['BAD_RE_ERR']
pa_shift = hdu[1].data['PA_SHIFT']
pa_shift_error = hdu[1].data['PA_SHIFT_ERR']
kine_twist = hdu[1].data['KINE_TWIST']
kine_twist_err = hdu[1].data['KINE_TWIST_ERR']
disturbed_kine = hdu[1].data['DISTURBED_KINE']
disturbed_kine_err = hdu[1].data['DISTURBED_KINE_ERR']
merging = hdu[1].data['MERGING']
merging_err = hdu[1].data['MERGING_ERR']
sym_OH = hdu[1].data['SYMMETRIC_OH']
sym_OH_err = hdu[1].data['SYMMETRIC_OH_ERR']
distorted_OH = hdu[1].data['DISTORTED_OH']
distorted_OH_err = hdu[1].data['DISTORTED_OH_ERR']
chaotic_OH = hdu[1].data['CHAOTIC_OH']
chaotic_OH_err = hdu[1].data['CHAOTIC_OH_ERR']
bad_OH = hdu[1].data['BAD_OH']
bad_OH_err = hdu[1].data['BAD_OH_ERR']
low_knots = hdu[1].data['LOW_KNOTS'] 
low_knots_err = hdu[1].data['LOW_KNOTS_ERR']
high_knots = hdu[1].data['HIGH_KNOTS']
high_knots_err = hdu[1].data['HIGH_KNOTS_ERR']
linear_OHgrad = hdu[1].data['LINEAR_OHGRAD']
linear_OHgrad_err = hdu[1].data['LINEAR_OHGRAD_ERR']
slope_change = hdu[1].data['SLOPE_CHANGE']
slope_change_err = hdu[1].data['SLOPE_CHANE_ERR']
irr_OHgrad = hdu[1].data['IRREGULAR_OHGRAD']
irr_OHgrad_err = hdu[1].data['IRREGULAR_OHGRAD_ERR']
bad_OHgrad = hdu[1].data['BAD_OHGRAD']
bad_OHgrad_err = hdu[1].data['BAD_OHGRAD_ERR']

for i in range (0, len(disturbed_kine)-1):
    if disturbed_kine == 1:
        print(plateifu)