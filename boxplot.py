"""
Thanks Josh Hemann for the example
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon


# Generate some data from five different probability distributions,
# each with different characteristics. We want to play with how an IID
# bootstrap resample of the data preserves the distributional
# properties of the original sample, and a boxplot is one visual tool
# to make this assessment
def make_boxplot(data, experiments, boxColors, methods, title):

	numLabels = len(experiments)

	fig = plt.figure(figsize=(10,6))
	fig.canvas.set_window_title('A Boxplot Example')
	ax1 = fig.add_subplot(111)
	plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
	
	bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='red', marker='+')
	
	# Add a horizontal grid to the plot, but make it very light in color
	# so we can use it for reading data values but not be distracting
	ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
	
	# Hide these grid behind plot objects
	ax1.set_axisbelow(True)
	ax1.set_title(title)
	ax1.set_xlabel('Datasets')
	ax1.set_ylabel('Euclidian distance')
	
	# Now fill the boxes with desired colors
	#boxColors = ['darkkhaki','royalblue']
	numBoxes = numLabels*len(methods)
	medians = range(numBoxes)
	for i in range(numBoxes):
	  box = bp['boxes'][i]
	  boxX = []
	  boxY = []
	  for j in range(len(experiments)):
	      boxX.append(box.get_xdata()[j])
	      boxY.append(box.get_ydata()[j])
	  boxCoords = zip(boxX,boxY)
	  # Alternate between Dark Khaki and Royal Blue
	  k = i % len(methods)
	  boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
	  ax1.add_patch(boxPolygon)
	  # Now draw the median lines back over what we just filled in
	  med = bp['medians'][i]
	  medianX = []
	  medianY = []
	  for j in range(2):
	      medianX.append(med.get_xdata()[j])
	      medianY.append(med.get_ydata()[j])
	      plt.plot(medianX, medianY, 'k')
	      medians[i] = medianY[0]
	  # Finally, overplot the sample averages, with horizontal alignment
	  # in the center of each box
	  plt.plot([np.average(med.get_xdata())], [np.average(data[i])], color='w', marker='*', markeredgecolor='k')
	
	# Set the axes ranges and axes labels
	ax1.set_xlim(0.5, numBoxes+0.5)
	top = min(max([max(d) for d in data]), 5)
	bottom = 0
	ax1.set_ylim(bottom, top)
	xtickNames = plt.setp(ax1, xticklabels=np.repeat(experiments, len(methods)))
	plt.setp(xtickNames, rotation=45, fontsize=8)
	
	# Due to the Y-axis scale being different across samples, it can be
	# hard to compare differences in medians across the samples. Add upper
	# X-axis tick labels with the sample medians to aid in comparison
	# (just use two decimal places of precision)
	pos = np.arange(numBoxes)+1
	upperLabels = [str(np.round(s, 3)) for s in medians]

	#weights = ['roman']*len(methods)
	weights = []
			
	for i in range(len(experiments)):
		min_median = min(medians[i*len(methods):(i+1)*len(methods)])
		for j in range(len(methods)):
			if medians[i*len(methods) + j] == min_median:
				weights.append('bold')
			else:
				weights.append('roman')
	#weights = ['bold', 'semibold'] * len(methods)
	wi = 0
	for tick,label in zip(range(numBoxes),ax1.get_xticklabels()):
	   k = tick % len(methods)
	   ax1.text(pos[tick], top-(top*0.05), upperLabels[tick], horizontalalignment='center', size='x-small', weight=weights[wi], color=boxColors[k])
	   wi += 1
	   #k += 1

	
	# Finally, add a basic legend
	legend_Y = 0.17
	legend_X = 0.8
	for i in range(len(methods)):
		plt.figtext(legend_X, legend_Y - i*0.035, methods[i], backgroundcolor=boxColors[i], color='black', weight='roman', size='x-small')

	plt.figtext(legend_X, legend_Y - 0.005 - 0.035*(len(methods)), '*', color='white', backgroundcolor='silver', weight='roman', size='medium')
	plt.figtext(legend_X+0.015, legend_Y - 0.035*(len(methods)), ' Average Value', color='black', weight='roman', size='x-small')
	
	plt.show()
	
if __name__ == "__main__":
	experiments = ["Test%d"%(i) for i in range(5)]
	methods = ["Method A", "Method B", "Method C", "Method D"]
	data = [[0.4, 0.3, 0.1*i] for i in range(len(experiments)*len(methods))]
	boxColors = ["red", "blue", "royalblue", "darkkhaki"]
	make_boxplot(data, experiments, boxColors, methods)	
