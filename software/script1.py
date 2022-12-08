# set color palette
color = bokeh.palettes.Category10_10

# initialize figure instance
p = bokeh.plotting.figure(
    width = 1200,
    height = 300,
    x_axis_label = 'time (min)',
    y_axis_label = 'area (µm²)'
)

# draw normalized fluorescence signal for each trial
for i, c, s in zip(bac_list, color, bac_list_string):
    p.circle(
        x=i['time (min)'],
        y=i['area (µm²)'],
        color=c,
        legend_label=s
    )

# set legend options and show plot    
p.add_layout(p.legend[0], 'right')
p.legend.click_policy = "hide"
bokeh.io.show(p)
