# create a color iterator
color = itertools.cycle(palette)

# initialize figure instance
r = bokeh.plotting.figure(
    width = 1200,
    height = 300,
    x_axis_label = 'growth time',
    y_axis_label = 'area (µm²)'
)

# draw normalized fluorescence signal for each trial
for i, c in zip(ge_unique, color):
    tempdf = df.loc[df['growth event']==i]
    r.line(
        x=tempdf['growth_time'],
        y=tempdf['area (µm²)'],
        color=c,
        legend_label=str(i)
    )

# set legend options and show plot    
r.add_layout(r.legend[0], 'right')
r.legend.click_policy = "hide"
bokeh.io.show(r)
