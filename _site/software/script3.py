colors = colorcet.b_glasbey_category10

s = bokeh.plotting.figure(title = 'linear MLE',
    x_axis_label="a0",
    y_axis_label="k",
    frame_height=400,
    frame_width=400,
)

for trial, bs_reps in reps.items():
    # Extract contour lines in D-k_off plane.
    x_line, y_line = bebi103.viz.contour_lines_from_samples(
        x=bs_reps[:, -3], y=bs_reps[:, -2], levels=[0.95]
    )

    # Plot the contour lines with fill
    for x, y in zip(x_line, y_line):
        s.line(x, y, line_width=2, color=colors[trial], legend_label=f"growth event {trial}")
        s.patch(x, y, fill_color=colors[trial], alpha=0.3)
        
for trial, bs_reps in reps.items():
    x = [df_mle.loc[trial, "a0"]] * 2
    y = np.percentile(bs_reps[:, -2], [2.5, 97.5])
    s.line(x, y, line_width=3, color=colors[trial])

    y = [df_mle.loc[trial, "k"]] * 2
    x = np.percentile(bs_reps[:, -3], [2.5, 97.5])
    s.line(x, y, line_width=3, color=colors[trial])     

s.add_layout(s.legend[0], 'right')
s.legend.click_policy = "hide"
bokeh.io.show(s)
