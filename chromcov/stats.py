def breadth_at_depth(base_depth, bins=10):
	# Take the cumulative sum in reverse order
	# (start at highest depth)
	hist, bin_edges = np.histogram(base_depth, bins=bins)
	cum_breadth_rev = np.cumsum(hist[::-1])

	# Reverse
	cum_breadth = cum_breadth_rev[::-1]
	cum_breadth_pct = cum_breadth / cum_breadth[0]

	return breadth_hist, breadth_pcts, bin_edges