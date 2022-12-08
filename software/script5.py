v = bokeh.plotting.figure(title = 'model comparison',
    x_axis_label="AIC_linear",
    y_axis_label="AIC_exponential",
    frame_height=400,
    frame_width=400                    
)

for i in ge_unique:
    v.circle(
        x=results.iloc[i,7],
        y=results.iloc[i,8],
        color=colors[i],
        #legend_label=str(i)
    )
    
v.line(x = np.linspace(120,720,1000),
       y = np.linspace(120,720,1000),
       color = 'gray', 
       legend_label = 'AIC_lin==AIC_exp'
)

# set legend options and show plot    
v.add_layout(v.legend[0], 'right')
v.legend.click_policy = "hide"
bokeh.io.show(v)
