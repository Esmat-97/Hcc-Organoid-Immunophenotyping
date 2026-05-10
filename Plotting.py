import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import gseapy as gp



# --- قراءة البيانات ---
df = pd.read_csv("GSE182593_FPKM_Matrix_Gene_Level.txt", sep="\t")
df.set_index("gene_id", inplace=True)

# --- اختيار الأعمدة الرقمية فقط ---
numeric_df = df.select_dtypes(include=["number"])

# --- تنظيف البيانات ---
# استبدال أي NaN بالقيمة 0
numeric_df = numeric_df.fillna(0)






# --- PCA ---
scaler = StandardScaler()
X_scaled = scaler.fit_transform(numeric_df.T)

pca = PCA(n_components=2)
pca_result = pca.fit_transform(X_scaled)

# --- رسم النتائج ---
samples = numeric_df.columns
plt.figure(figsize=(8,6))
plt.scatter(pca_result[:,0], pca_result[:,1], alpha=0.7)
for i, sample in enumerate(samples):
    plt.text(pca_result[i,0], pca_result[i,1], sample, fontsize=9)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA Projection of PDO Samples")
plt.savefig("PCA_PDOs.png", dpi=300) 
plt.show()




# --- 6. Heatmap لأعلى الجينات المتغيرة ---
top_genes = numeric_df.var(axis=1).sort_values(ascending=False).head(50).index
plt.figure(figsize=(12,8))
sns.heatmap(numeric_df.loc[top_genes], cmap="viridis")
plt.title("Top 50 Variable Genes in PDOs")
plt.savefig("Top50_Variable_Genes.png", dpi=300)
plt.show()







# --- 4. استخراج Biomarkers أساسية ---
biomarkers = ["IFNG", "GZMB", "PRF1"]
subset = numeric_df.loc[numeric_df.index.intersection(biomarkers)]

plt.figure(figsize=(8,6))
sns.heatmap(subset, cmap="coolwarm", annot=True)
plt.title("Immune Biomarkers (IFNG, GZMB, PRF1)")
plt.savefig("Immune_Biomarkers.png", dpi=300)
plt.show()







immune_genes = ["IFNG","CXCL9","CXCL10","PRF1","GZMB","NKG7",
                "CD8A","PDCD1","LAG3","TIGIT","HAVCR2"]

immune_subset = numeric_df.loc[numeric_df.index.intersection(immune_genes)]

plt.figure(figsize=(12,6))
sns.heatmap(immune_subset, cmap="coolwarm", annot=False)
plt.title("Immune-focused Heatmap (Hot vs Cold PDOs)")
plt.savefig("Immune_Focused.png", dpi=300)
plt.show()










# --- Immune Score جديد (5 جينات) ---
responder_genes = ["IFNG","GZMB","PRF1","CXCL9","CXCL10"]
immune_score = numeric_df.loc[numeric_df.index.intersection(responder_genes)].sum()
plt.figure(figsize=(10,6))
sns.barplot(x=immune_score.index, y=immune_score.values, palette="coolwarm")
plt.xticks(rotation=90)
plt.ylabel("immune score")
plt.title("immune score per Sample")
plt.savefig("immune score.png", dpi=300)
plt.show()
print("\n=== Immune Score per Sample ===")
print(immune_score)






exhaustion_genes = ["PDCD1","LAG3","TIGIT","HAVCR2","CTLA4"]
exhaustion_score = numeric_df.loc[numeric_df.index.intersection(exhaustion_genes)].sum()
plt.figure(figsize=(10,6))
sns.barplot(x=exhaustion_score.index, y=exhaustion_score.values, palette="coolwarm")
plt.xticks(rotation=90)
plt.ylabel("exhaustion score")
plt.title("exhaustion score per Sample")
plt.savefig("exhaustion score.png", dpi=300)
plt.show()
print("\n=== Exhaustion Score per Sample ===")
print(exhaustion_score)






cytotoxic_genes = ["IFNG","GZMB","PRF1","NKG7","GNLY"]
cytotoxic_score = numeric_df.loc[numeric_df.index.intersection(cytotoxic_genes)].sum()
plt.figure(figsize=(10,6))
sns.barplot(x=cytotoxic_score.index, y=cytotoxic_score.values, palette="coolwarm")
plt.xticks(rotation=90)
plt.ylabel("cytotoxic score ")
plt.title("cytotoxic score  per Sample")
plt.savefig("cytotoxic score.png", dpi=300)
plt.show()
print("\n=== Cytotoxic Score per Sample ===")
print(cytotoxic_score)




antigen_genes = ["HLA-A","HLA-B","B2M","TAP1","TAP2","NLRC5"]
antigen_score = numeric_df.loc[numeric_df.index.intersection(antigen_genes)].sum()
plt.figure(figsize=(10,6))
sns.barplot(x=antigen_score.index, y=antigen_score.values, palette="coolwarm")
plt.xticks(rotation=90)
plt.ylabel("antigen_score ")
plt.title("antigen_score  per Sample")
plt.savefig("antigen_score.png", dpi=300)
plt.show()
print("\n=== Antigen Presentation Score per Sample ===")
print(antigen_score)





chemokines = ["CCL5","CXCR3","CCR5"]
chemokine_score = numeric_df.loc[numeric_df.index.intersection(chemokines)].sum()
plt.figure(figsize=(10,6))
sns.barplot(x=chemokine_score.index, y=chemokine_score.values, palette="coolwarm")
plt.xticks(rotation=90)
plt.ylabel("chemokine_score ")
plt.title("chemokine_score  per Sample")
plt.savefig("chemokine_score.png", dpi=300)
plt.show()
print("\n=== Chemokine Recruitment Score per Sample ===")
print(chemokine_score)





nk_genes = ["NKG7","GNLY","NCR1","KLRD1","FCGR3A","MICA","MICB"]
nk_score = numeric_df.loc[numeric_df.index.intersection(nk_genes)].sum()
plt.figure(figsize=(10,6))
sns.barplot(x=nk_score.index, y=nk_score.values, palette="coolwarm")
plt.xticks(rotation=90)
plt.ylabel("nk_score ")
plt.title("nk_score  per Sample")
plt.savefig("nk_score.png", dpi=300)
plt.show()
print("\n=== NK Activation Score per Sample ===")
print(nk_score)




stromal_genes = ["COL1A1","ACTA2","FAP","TGFB1","VIM"]
stromal_score = numeric_df.loc[numeric_df.index.intersection(stromal_genes)].sum()
plt.figure(figsize=(10,6))
sns.barplot(x=stromal_score.index, y=stromal_score.values, palette="coolwarm")
plt.xticks(rotation=90)
plt.ylabel("stromal_score ")
plt.title("stromal_score  per Sample")
plt.savefig("stromal_score.png", dpi=300)
plt.show()
print("\n=== Stromal/Fibrotic Score per Sample ===")
print(stromal_score)




emt_genes = ["EPCAM","CD44","PROM1","SOX9","VIM","ZEB1"]
emt_score = numeric_df.loc[numeric_df.index.intersection(emt_genes)].sum()
plt.figure(figsize=(10,6))
sns.barplot(x=emt_score.index, y=emt_score.values, palette="coolwarm")
plt.xticks(rotation=90)
plt.ylabel("emt_score")
plt.title("emt_score  per Sample")
plt.savefig("emt_score.png", dpi=300)
plt.show()
print("\n=== EMT/Stemness Score per Sample ===")
print(emt_score)







# --- تعريف المسارات المناعية ---
pathways = {
    "IFN_gamma_response": ["IFNG","CXCL9","CXCL10"],
    "Cytotoxic_Tcell": ["CD8A","GZMB","PRF1","NKG7"],
    "NK_activation": ["NKG7","KLRD1","GNLY"],
    "Antigen_presentation": ["HLA-A","HLA-B","HLA-C","B2M"],
    "Exhaustion": ["PDCD1","LAG3","TIGIT","HAVCR2"]
}

# --- حساب الـ scores لكل pathway ---
scores = {}
for pw, genes in pathways.items():
    # نجمع التعبير الجيني للجينات المختارة
    scores[pw] = numeric_df.loc[numeric_df.index.intersection(genes)].sum()

# تحويل النتائج إلى DataFrame مرتب
scores_df = pd.DataFrame(scores)

# --- رسم Heatmap ---
plt.figure(figsize=(12,6))
sns.heatmap(scores_df.T, cmap="coolwarm", annot=False)
plt.title("Immune Pathway Scores (PDOs)")
plt.xlabel("Samples")
plt.ylabel("Immune Pathways")
plt.savefig("Immune_Pathway_Scores.png", dpi=300) 
plt.show()








import matplotlib.pyplot as plt

# افترض إن عندك Series cytotoxic_score و exhaustion_score
plt.figure(figsize=(8,6))
plt.scatter(exhaustion_score, cytotoxic_score, alpha=0.7)

for sample in cytotoxic_score.index:
    plt.text(exhaustion_score[sample], cytotoxic_score[sample], sample, fontsize=8)

plt.xlabel("Exhaustion Score")
plt.ylabel("Cytotoxic Score")
plt.title("Immune Classification Scatter Plot")

plt.axhline(y=0.1, color='grey', linestyle='--')  # خط إرشادي
plt.axvline(x=0.5, color='grey', linestyle='--')  # خط إرشادي
plt.savefig("Immune_Classification.png", dpi=300)
plt.show()






import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore

# --- نجمع كل الـ scores في DataFrame واحد ---
all_scores = pd.DataFrame({
    "Immune": immune_score,
    "Exhaustion": exhaustion_score,
    "Cytotoxic": cytotoxic_score,
    "Chemokine": chemokine_score,
    "NK": nk_score,
    "EMT": emt_score,
    "Stromal": stromal_score,
    "Antigen": antigen_score
})

print("\n=== Integrated Scores Table ===")
print(all_scores)

# --- z-score normalization ---
all_scores_z = all_scores.apply(zscore)

# --- Hierarchical clustering heatmap ---
plt.figure(figsize=(12,8))
sns.clustermap(all_scores_z,
               cmap="coolwarm",
               method="ward",      # طريقة clustering
               metric="euclidean", # المسافة
               standard_scale=None,
               figsize=(12,8),
               row_cluster=True,
               col_cluster=True)

plt.title("Integrated Phenotype Heatmap (PDOs)", pad=80)
plt.savefig("Integrated_Phenotype_Heatmap.png", dpi=300)
plt.show()
