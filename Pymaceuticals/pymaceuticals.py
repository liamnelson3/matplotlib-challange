# %% [markdown]
# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# - Add your analysis here.
#  

# %%
# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)



# %%
# Combine the data into a single DataFrame
merged_mouse_df = pd.merge(study_results, mouse_metadata, on = "Mouse ID", how = "left")

# Display the data table for preview
merged_mouse_df.head()

# %%
# Checking the number of mice.
number_of_mice = merged_mouse_df["Mouse ID"].nunique()
number_of_mice
#len(merged_mouse_df["Mouse ID"].unique()) works to

# %%
# Our data should be uniquely identified by Mouse ID and Timepoint
# Get the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mice = merged_mouse_df[merged_mouse_df.duplicated(["Mouse ID", "Timepoint"])]
duplicate_mice = duplicate_mice["Mouse ID"].unique()
print(duplicate_mice)

# %%
# Optional: Get all the data for the duplicate mouse ID. 
duplicate_df=merged_mouse_df.loc[(merged_mouse_df["Mouse ID"] == "g989")]
duplicate_df

# %%
# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_mice_df = merged_mouse_df.loc[(merged_mouse_df["Mouse ID"] != "g989")]
clean_mice_df

# %%
# Checking the number of mice in the clean DataFrame.
clean_mice_df["Mouse ID"].nunique()

# %% [markdown]
# ## Summary Statistics

# %%
# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.
drug_regimen_mean = clean_mice_df.groupby("Drug Regimen")["Tumor Volume (mm3)"].mean()
drug_regimen_median = clean_mice_df.groupby("Drug Regimen")["Tumor Volume (mm3)"].median()
drug_regimen_variance = clean_mice_df.groupby("Drug Regimen")["Tumor Volume (mm3)"].var()
drug_regimen_std = clean_mice_df.groupby("Drug Regimen")["Tumor Volume (mm3)"].std()
drug_regimen_sem = clean_mice_df.groupby("Drug Regimen")["Tumor Volume (mm3)"].sem()

drug_regimen_df = pd.merge(pd.merge(pd.merge(pd.merge(drug_regimen_mean, drug_regimen_median, on = "Drug Regimen"), 
                                                      drug_regimen_variance, on = "Drug Regimen"), 
                                                      drug_regimen_std, on = "Drug Regimen"), 
                                                      drug_regimen_sem, on = "Drug Regimen")

drug_regimen_df.columns.values[0] = "Mean Tumor Volume"
drug_regimen_df.columns.values[1] = "Median Tumor Volume"
drug_regimen_df.columns.values[2] = "Tumor Volume Variance"
drug_regimen_df.columns.values[3] = "Tumor Volume Std. Dev."
drug_regimen_df.columns.values[4] = "Tumor Volume Std. Err."

drug_regimen_df

# %%
# A more advanced method to generate a summary statistics table of mean, median, variance, standard deviation,
# and SEM of the tumor volume for each regimen (only one method is required in the solution)

# Using the aggregation method, produce the same summary statistics in a single line
aggregated_mice_df = clean_mice_df.groupby("Drug Regimen")["Tumor Volume (mm3)"].aggregate(['mean', 'median', 'var', 'std', 'sem'])
aggregated_mice_df

# %% [markdown]
# ## Bar and Pie Charts

# %%
#setup
mice_drug_groups = clean_mice_df.groupby("Drug Regimen")
mice_drug_counts = mice_drug_groups[["Timepoint"]].count().sort_values("Timepoint", ascending = False)

# %%
# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using Pandas.
mice_drug_counts.plot(kind = "bar", ylabel = "# of Observed Mouse Timestamps")
plt.show

# %%
# Generate a bar plot showing the total number of rows (Mouse ID/Timepoints) for each drug regimen using pyplot.
mice_drug_counts = mice_drug_counts.reset_index()

drug_regimen_bars = mice_drug_counts["Drug Regimen"]
mouse_timestamps_bars = mice_drug_counts["Timepoint"]

plt.bar(drug_regimen_bars, mouse_timestamps_bars, color="royalblue",  width=0.5, align="center")
plt.xlabel("Drug Regimen")
plt.ylabel("# of Observed Mouse Timestamps")
plt.xticks(rotation = "vertical")
plt.show()

# %%
#setup
mice_gender_ratio = clean_mice_df.groupby("Sex")[["Mouse ID"]].count()
mice_gender_ratio

# %%
# Generate a pie plot showing the distribution of female versus male mice using Pandas
mice_gender_ratio.plot(kind = "pie", subplots = True, autopct='%1.0f%%')

# %%
# Generate a pie plot showing the distribution of female versus male mice using pyplot
mice_gender_ratio = mice_gender_ratio.reset_index()
sizes = mice_gender_ratio["Mouse ID"]
labels = mice_gender_ratio["Sex"]

plt.pie(sizes, labels = labels, autopct='%1.0f%%')

# %% [markdown]
# ## Quartiles, Outliers and Boxplots

# %%
# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin
options = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']

#filters rows where Drug Regimen equals one of the items in the list options
selective_drug_df = clean_mice_df[clean_mice_df['Drug Regimen'].isin(options)]
selective_drug_df.sort_values("Mouse ID")

# %%
# Start by getting the last (greatest) timepoint for each mouse
selective_drug_groups = selective_drug_df.groupby( "Mouse ID")

max_timepoints = selective_drug_groups[["Timepoint"]].max()
max_timepoints.reset_index()
# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
max_timepoints_merged = pd.merge(merged_mouse_df, max_timepoints, on=("Mouse ID", "Timepoint"), how = "right")
max_timepoints_merged.sort_values("Mouse ID")

# %%
# Put treatments into a list for for loop (and later for plot labels)
drug_treatment_names = options

# Create empty list to fill with tumor vol data (for plotting)
tumor_vol_list = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
for i in drug_treatment_names:
    # Locate the rows which contain mice on each drug and get the tumor volumes
    tumor_vol_iter = max_timepoints_merged.loc[max_timepoints_merged["Drug Regimen"] == i, "Tumor Volume (mm3)"]

    # add subset
    tumor_vol_list.append(tumor_vol_iter)

    # Determine outliers using upper and lower bounds
    quartiles = tumor_vol_iter.quantile([.25,.5,.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq-lowerq

    lower_bound = lowerq - (1.5*iqr)
    upper_bound = upperq + (1.5*iqr)
    outliers = tumor_vol_iter.loc[(tumor_vol_iter < lower_bound) | 
                                  (tumor_vol_iter > upper_bound)]
    print(f"For {i}, the outliers are: {outliers}")
#dont worry about the number before the 36.321346. thats just the index number for the outlier

# %%
# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
fig1, ax1 = plt.subplots()
ax1.set_ylabel('Final Tumor Volume (mm3)')
ax1.boxplot(tumor_vol_list, labels = drug_treatment_names)
plt.show()

# %% [markdown]
# ## Line and Scatter Plots

# %%
# Generate a line plot of tumor volume vs. time point for a single mouse treated with Capomulin
capomulin_regimen = clean_mice_df.loc[clean_mice_df["Drug Regimen"] == "Capomulin"]
mouse_i509 = capomulin_regimen.loc[capomulin_regimen["Mouse ID"] == "l509"]
x_axis = mouse_i509["Timepoint"]
y_axis = mouse_i509["Tumor Volume (mm3)"]
plt.xlabel("Capomulin treatment of Mouse l509")
plt.xlabel("Tumor Volume (mm3)")
plt.ylabel("Timepoint (Day)")
plt.plot(x_axis, y_axis)

# %%
# Generate a scatter plot of mouse weight vs. the average observed tumor volume for the entire Capomulin regimen
capomulin_regimen_average = capomulin_regimen.groupby(["Mouse ID"]).mean()

x_weight = capomulin_regimen_average["Weight (g)"]
y_tumor = capomulin_regimen_average["Tumor Volume (mm3)"]

plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")
plt.scatter(x_weight, y_tumor)

# %% [markdown]
# ## Correlation and Regression

# %%
# Calculate the correlation coefficient and a linear regression model 
# for mouse weight and average observed tumor volume for the entire Capomulin regimen

# Correlation Coefficient
pearson_coef = st.pearsonr(x_weight, y_tumor)
print(f"The correlation is {round(pearson_coef[0],2)}")

# Linear Regression
slope, intercept, rvalue, pvalue, stderr = st.linregress(x_weight, y_tumor)
regress_values = x_weight * slope + intercept

plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")

plt.scatter(x_weight, y_tumor)
plt.plot(x_weight,regress_values,"r-")
plt.show()

# %%



