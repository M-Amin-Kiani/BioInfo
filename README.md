
# 🧬 BioInfo

## 📌 درباره‌ی پروژه

پروژه‌ی **BioInfo** شامل اسکریپت‌ها و نوت‌بوک‌هایی برای تحلیل داده‌های بیوانفورماتیکی است. این تحلیل‌ها با استفاده از زبان‌های پایتون و R انجام می‌شوند.

---

## 🛠️ پیش‌نیازها

- زبان‌های برنامه‌نویسی:
  - Python
  - R

- کتابخانه‌های پایتون:
  - pandas
  - numpy
  - scipy
  - biopython

- کتابخانه‌های R:
  - tidyverse
  - Bioconductor

### نصب کتابخانه‌های پایتون:

```bash
pip install pandas numpy scipy biopython
```

### نصب کتابخانه‌های R:

```R
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

---

## 📂 ساختار پروژه

```
BioInfo/
├── research/
│   ├── analysis1.ipynb
│   ├── analysis2.R
│   └── data/
│       ├── dataset1.csv
│       └── dataset2.fasta
└── README.md
```

---

## 🚀 نحوه‌ی اجرا

1. کلون کردن ریپازیتوری:

```bash
git clone https://github.com/M-Amin-Kiani/BioInfo.git
cd BioInfo
```

2. نصب پیش‌نیازها (مطابق بالا)

3. اجرای نوت‌بوک یا اسکریپت:

- اجرای نوت‌بوک Jupyter:

```bash
jupyter notebook research/analysis1.ipynb
```

- اجرای اسکریپت R:

```bash
Rscript research/analysis2.R
```

---

## ✍️ نویسنده

- محمد امین کیانی  
- GitHub: [M-Amin-Kiani](https://github.com/M-Amin-Kiani)

---

## 📄 مجوز

این پروژه تحت مجوز MIT منتشر شده است.

---

## 📬 پشتیبانی

در صورت وجود سوال یا مشکل، لطفاً issue جدیدی در گیت‌هاب باز کنید.
