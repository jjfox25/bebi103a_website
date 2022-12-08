results['w_linear'] = np.exp(-results['AIC_lin']/2) / (np.exp(-results['AIC_lin']/2) + np.exp(-results['AIC_exp']/2))
w_lin_striphist = iqplot.stripbox(data=results, q="w_linear", title="AIC Weight linear", spread="jitter", q_axis='y')

bokeh.io.show(w_lin_striphist)
