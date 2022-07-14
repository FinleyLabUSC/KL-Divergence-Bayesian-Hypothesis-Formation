import matplotlib.pyplot as plt
import numpy as np

fontsize=14

#create supplementary figure S3A, B, and C

#re-plot Jiao 2013 Figure 6 - data from plot digitizer
labelvector = ['Cish','Socs2','Socs3','Socs4','Socs5','Socs6','Socs7']
controlexpression = [0.9613259668508287, 0.9613259668508287, 0.9613259668508287, 0.9613259668508287, 0.9613259668508287, 0.9613259668508287, 0.9613259668508287]
pregnantexpression = y = [2.61878453038674, 2.187845303867403, 0.6464088397790055, 0.6961325966850829, 0.8287292817679558, 0.8121546961325966, 0.6795580110497237]


x = np.arange(len(labelvector))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, controlexpression, width, color='silver')
rects2 = ax.bar(x + width/2, pregnantexpression, width, color='dimgrey')
rects3 = ax.bar(x[0:2] - width/2, controlexpression[0:2], width, label='Control Group', color='lightcoral')
rects4 = ax.bar(x[0:2] + width/2, pregnantexpression[0:2], width, label='Pregnant Group', color='firebrick')



# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Relative Expression (fold)',fontsize=fontsize)
#ax.set_title('Jiao 2013, Figure 6c')
ax.set_xticks(x)
ax.set_xticklabels(labelvector, fontsize=fontsize)
ax.legend(fontsize=fontsize)


fig.tight_layout()

plt.show()
plt.savefig('literatureevidence/Jiao2013Figure6.png', dpi=300)

#re-plot Rieck 2009 Figure 3 - data from plot digitizer
day = [0, 10.5, 14.5, 18.5]
cishexpression = [1.0185185185185184, 4.0740740740740735, 2.7777777777777777, 2.8703703703703702]
socs2expression = [0.9818181818181818, 1.9636363636363636, 2.1818181818181817, 2.018181818181818]
x = np.arange(len(day)) # the label locations

plt.subplot(2, 1, 1)
plt.bar(x, cishexpression, color='firebrick')
plt.xticks(x, labels=day)
plt.ylim([0,5])
plt.ylabel('Cish Expression',fontsize=fontsize)
plt.subplot(2, 1, 2)
plt.bar(x, socs2expression, color='firebrick')
plt.xticks(x, labels=day)
plt.ylabel('Socs2 Expression',fontsize=fontsize)
plt.ylim([0, 3])
plt.yticks([0, 1, 2, 3])
plt.xlabel('Day of Pregnancy',fontsize=fontsize)
plt.show()
plt.savefig('literatureevidence/Rieck2009Figure3.png', dpi=300)


#re-plot Chong 2001 Figure 2 - data from plot digitizer
time = [0, 1, 2, 4, 24]
CISH_TNFA = [1, 1, 6.237288135593221, 4.88135593220339, 3.6610169491525424, 4.338983050847458]
CISH_IL1B = [1, 3.2542372881355934, 8.542372881355933, 3.7966101694915255, 1.4915254237288136]
CISH_IFNG = [1, 4.610169491525424, 5.830508474576271, 4.067796610169491, 1.8983050847457628]

SOCS1IFNG = [1, 1.7959183673469385, 5.387755102040816, 8.653061224489795, 8.326530612244897]
SOCS1IL1B = [1, 0.9795918367346939, 1.4693877551020407, 0.7346938775510203, 0.9795918367346939]
SOCS1TNFA = [1, 0.9795918367346939, 1.4693877551020407, 0.9795918367346939, 0.9795918367346939]

SOCS2IL1B = [1, 1.3333333333333333, 17.333333333333332, 11, 3.333333333333333]
SOCS2TNFA = [1, 8, 12.666666666666666, 3, 1.6666666666666665]
SOCS2IFNG = [1, 10.545454545454545, 7.272727272727273, 6.545454545454546, 1.0909090909090908]

plt.figure(figsize=(13, 4))
plt.subplots_adjust(bottom=0.15, left=0.05, right=0.90)
plt.subplot(1,3,1)
plt.plot(time, CISH_IFNG, 'o-', color='black', label='CISH', markersize=8)
plt.xticks([0, 1, 2, 4, 24], labels=time)
plt.xlabel('Time, Hours', fontsize=fontsize)
plt.ylabel('Relative Expression', fontsize=fontsize)
plt.legend(loc='center right', fontsize=fontsize)

plt.subplots_adjust(bottom=0.15, top=0.99, left=0.05, right=0.95)
plt.subplot(1,3,2)
plt.plot(time, SOCS1IFNG, 'o-', color='black', label='SOCS1', markersize=8)
plt.xticks([0, 1, 2, 4, 24], labels=time)
plt.xlabel('Time, Hours', fontsize=fontsize)
#plt.ylabel('SOCS1 Expression')
plt.legend(loc='center right',fontsize=fontsize)

plt.subplots_adjust(bottom=0.15, top=0.99, left=0.05, right=0.95)
plt.subplot(1,3,3)
plt.plot(time, SOCS2IFNG, 'o-', color='black', label='SOCS2', markersize=8)
plt.xticks([0, 1, 2, 4, 24], labels=time)
#plt.xlabel('Time, Hours')
#plt.ylabel('SOCS2 Expression')
plt.legend(loc='center right',fontsize=fontsize)

plt.show()
plt.savefig('literatureevidence/Chong2001Figure2.png', dpi=300)
