results['w_exp'] = np.exp(-results['AIC_exp']/2) / (np.exp(-results['AIC_exp']/2) + np.exp(-results['AIC_lin']/2))

w_exp_striphist = iqplot.stripbox(data=results, q="w_exp", title="AIC Weight exponential", spread="jitter", q_axis='y')

bokeh.io.show(w_exp_striphist)
