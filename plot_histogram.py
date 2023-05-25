import pandas as pd
import matplotlib.pyplot as plt

# Read the data from a csv file
data = pd.read_csv('output.csv')

# Plot the data
plt.bar(data['bin'], data['amount'])

# Set plot labels and title
plt.xlabel('Bin')
plt.ylabel('Amount')
plt.title('Bins Histogram')

# Show the plot
plt.show()