def fit_gp_george_2d(spec, phases, waves, printlog=True,plot=True):
    import george
    from george import kernels
    import numpy as np
    import matplotlib.pyplot as plt

    # scale = 1.e13
    scale = 1.
    inflate_var = 10.
    y = spec[:,2]*scale
    yerr = inflate_var*scale/np.sqrt(spec[:,3])
    x = spec[:,0:2]
    
    idx = ~(np.isnan(y) | np.isnan(x[:,0]) | np.isnan(x[:,1]) | np.isnan(yerr))
    y = y[idx]
    yerr = yerr[idx]
    x = x[idx,:]
    # print(x)
    # print(y)
    # print(yerr)

    amp = np.var(y)
    tlengthscale = 10.
    wlengthscale = 10.*np.average(np.diff(waves))
    print(tlengthscale, wlengthscale)
#     kernel = amp * kernels.Matern32Kernel(ndim=2,metric=[tlengthscale,wlengthscale])
    kernel = amp * kernels.ExpSquaredKernel(ndim=2,metric=[tlengthscale,wlengthscale])
#     kernel.freeze_parameter("k2:metric:log_M_1_1")
    gp = george.GP(kernel,mean=np.median(y))
    gp.compute(x, yerr)

    print(gp.get_parameter_names())
    print(gp.get_parameter_vector())
    if printlog:
        print("Initial ln-likelihood: {0:.2f}".format(gp.log_likelihood(y)))

    def neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.log_likelihood(y)

    def grad_neg_ln_like(p):
        gp.set_parameter_vector(p)
        return -gp.grad_log_likelihood(y)

    from scipy.optimize import minimize 
    result = minimize(neg_ln_like, gp.get_parameter_vector(), jac=grad_neg_ln_like)
    print(result)

    gp.set_parameter_vector(result.x)
    if printlog:
        print("Final ln-likelihood: {0:.2f}".format(gp.log_likelihood(y)))

#     meanlc = []
#     lcvar = []
#     samples = []
#     peakmjds = []
#     filters = []
#     sample_funcs = []

    if plot:
        pind = int(len(phases)/2.)
        phasearr = [phases[0],phases[pind],phases[-1]]
        for p in phasearr:
            wave = np.linspace(waves[0],waves[-1],200)
            xa = np.vstack([[p]*200,wave]).T
    #         print(xa)
            pred_f, pred_var_f = gp.predict(y, xa, return_var=True)

    #         print(pred_f)
            xphase = x[:,0]
            idx = np.array(xphase==p)
    #         print(x[idx,1])
    #         print(y[idx])
            # print(x[idx,1])
            # print(y[idx])
            # print(yerr[idx])
            # plt.errorbar(x[idx,1],y[idx], yerr=yerr[idx], fmt='o')
            plt.errorbar(x[idx,1],y[idx], yerr=yerr[idx], fmt='o')
            plt.plot(wave,pred_f,color='k')
            plt.fill_between(wave,pred_f-np.sqrt(pred_var_f),pred_f+np.sqrt(pred_var_f),alpha=0.5,color='k')
            plt.show()

        wind1 = int(len(waves)/4.)
        wind2 = int(len(waves)/2.)
        wind3 = int(3*len(waves)/4.)
        wavearr = [waves[wind1],waves[wind2],waves[wind3]]
        for wave in wavearr:
            phasearr = np.linspace(-10,40,100)
            xa = np.vstack([phasearr,[wave]*100]).T
    #         print(xa)
            pred_f, pred_var_f = gp.predict(y, xa, return_var=True)

    #         print(pred_f)
            xwave = x[:,1]
            idx = np.array(xwave==wave)
            # print(x[idx,0])
            # print(y[idx])
            # print(yerr[idx])
            plt.errorbar(x[idx,0],y[idx], yerr=yerr[idx], fmt='o')
            plt.plot(phasearr,pred_f,color='k')
            plt.fill_between(phasearr,pred_f-np.sqrt(pred_var_f),pred_f+np.sqrt(pred_var_f),alpha=0.5,color='k')
            plt.show()
            
        #2d plot

        
    return gp, x, y
#         if nrandomsamples>0:
#             samples_f = gp.sample_conditional(y, xa, size=nrandomsamples)
#             samples.append(samples_f)
#             sample_func_f = []
#             for j in range(nrandomsamples):
#                 finterp = interp1d(xa0,samples_f[j],fill_value="extrapolate")
#                 sample_func_f.append(finterp)
#             sample_funcs.append(sample_func_f)
#         else:
#             sample_funcs = None
        
#         ##estimate peak in each band
#         meanf = interp1d(xa0,pred_f,fill_value="extrapolate")
#         varf = interp1d(xa0,pred_var_f,fill_value="extrapolate")
#         peakres = minimize(meanf,x[np.argmin(y[idx]),0])
#         pmjd = peakres.x[0]
#         if pmjd < x[idx,0].min() or pmjd > x[idx,0].max():
#             if printlog:
#                 print("No data before or after the peak for band", f)
#             pmjd = -99.
#         peakmjds.append(pmjd)
        
#         filters.append(f)
#         meanlc.append(meanf)
#         lcvar.append(varf)

#         if plot:
#             fig.add_subplot(len(sn['Filter'].unique())//6+1,6,i+1)
#             for j in range(nrandomsamples):
#                 plt.plot(xa[:,0],samples_f[j,:],'-',lw=0.5,color='g',alpha=0.5)
#             plt.fill_between(xa[:,0], pred_f - np.sqrt(pred_var_f), pred_f + np.sqrt(pred_var_f),
#                             color="k", alpha=0.2)
#             plt.plot(xa[:,0], pred_f, "k", lw=1.5, alpha=0.5)
#             plt.errorbar(x[idx,0], y[idx], yerr=yerr[idx], fmt=".k", capsize=0)
#             plt.ylim(plt.ylim()[::-1])
#             plt.title(f)
#             i += 1

#     gppars = gp.get_parameter_dict(include_frozen=True)
#     return GPResults(meanlc,lcvar,peakmjds,filters,
#                      gppars=gppars,realizations=sample_funcs,
#                      chisq=gp.log_likelihood(y))