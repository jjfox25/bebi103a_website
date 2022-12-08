u = bokeh.plotting.figure(title = 'exponential MLE',
    x_axis_label="a0",
    y_axis_label="k",
    frame_height=400,
    frame_width=400,
)

for trial, bs_reps in reps_exp.items():
    # Extract contour lines in D-k_off plane.
    x_line, y_line = bebi103.viz.contour_lines_from_samples(
        x=bs_reps[:, -3], y=bs_reps[:, -2], levels=[0.95]
    )

    # Plot the contour lines with fill
    for x, y in zip(x_line, y_line):
        u.line(x, y, line_width=2, color=colors[trial], legend_label=f"growth event {trial}")
        u.patch(x, y, fill_color=colors[trial], alpha=0.3)
        
for trial, bs_reps in reps_exp.items():
    x = [df_mle_exp.loc[trial, "a0"]] * 2
    y = np.percentile(bs_reps[:, -2], [2.5, 97.5])
    u.line(x, y, line_width=3, color=colors[trial])

    y = [df_mle_exp.loc[trial, "k"]] * 2
    x = np.percentile(bs_reps[:, -3], [2.5, 97.5])
    u.line(x, y, line_width=3, color=colors[trial])     

u.add_layout(u.legend[0], 'right')
u.legend.click_policy = "hide"
bokeh.io.show(u)
