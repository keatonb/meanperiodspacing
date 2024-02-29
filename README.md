# mean period spacing
Statistical tests for mean spacings among a set of values.

If you use this code in your work, I would appreciate an acknowledgement.

I describe how these tests work in a blog post at https://keatonb.github.io/archivers/mps

Usage example:
```python
# Values to test for mean spacings
values = np.array([640.533, 678.433, 597.554, 485.527, 604.642, 557.649, 749.27, 865.97]) 

# Determine what spacings to sample in the range of interest
maxperiod = 60
minperiod = 15
periodsample = getperiodsampling(values, minperiod, maxperiod, oversample_factor=20)

# Run all three mean-period-spacing tests
iv = IV(values, periodsample)
ks = KS(values, periodsample)
ft = FT(values, periodsample)

# Plot the results
fig, axs = plt.subplots(3,1)

axs[0].plot(periodsample, iv)
axs[0].set_title("Inverse Variance")
axs[0].set_ylabel("IV")
axs[0].set_ylim(0)

axs[1].plot(periodsample, ks)
axs[1].set_title("K-S Test")
axs[1].set_ylabel("log(Q)")

axs[2].plot(periodsample, ft)
axs[2].set_title("Fourier Transform")
axs[2].set_ylabel("power")
axs[2].set_ylim(0)

for ax in axs:
    ax.set_xlim(periodsample[0],periodsample[-1])
    ax.set_xlabel("period")

plt.tight_layout()
plt.show()
```
Here is the resulting plot, with all three tests supporting a mean period spacing around 38 seconds.

![mpstests](https://github.com/keatonb/meanperiodspacing/assets/6413923/4eae78b8-6cf6-431a-9fae-17d16a841733)
