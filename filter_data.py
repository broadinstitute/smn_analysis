import pandas as pd
import matplotlib.pyplot as plt

def pltcolor(sma_status):
    sma_colors = []
    for sample in sma_status:
        if sample=='Unknown':
            sma_colors.append('grey')
        elif sample=='Not Affected':
            sma_colors.append('blue')
        else:
            sma_colors.append('red')
    return sma_colors

df = pd.read_csv("example.tsv", sep='\t')
smn_c_count = df['c_count']
smn_total = df['total_count']
sma_status = df['SMA_status']
x = list(smn_c_count)
y = list(smn_total)
sma_colors = pltcolor(sma_status)
plt.xlim(-5, 700)
plt.ylim(-5, 700)
plt.scatter(x=x,y=y,c = sma_colors, marker='+', alpha=0.5)
plt.xlabel('SMN Reads with C')
plt.ylabel('SMN c.840 position Total Reads')
plt.title('SMN Reads with C  vs. Total Reads at c.840 Position')
plt.savefig('C_vs_total.png', dpi = 300)