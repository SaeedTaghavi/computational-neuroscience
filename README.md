# کتاب علوم اعصاب محاسباتی

> An online book on computational neuroscience, written in Persian — from excitable membranes to whole-brain modeling, with a physics and dynamical-systems perspective.

یک کتاب آنلاین  به زبان فارسی دربارهٔ علوم اعصاب محاسباتی — از سازوکارهای پایهٔ نورون تا مدل‌سازی شبکه‌ها و مغز کامل، با نگاهی از جنس فیزیک و ریاضی.

🌐 **مطالعهٔ آنلاین:** [computational-neuroscience.ir](https://computational-neuroscience.ir)

## دربارهٔ کتاب

این کتاب رفتار نورون‌ها و شبکه‌های عصبی را با زبان ریاضی و مدل‌سازی توصیف می‌کند و گام‌به‌گام پیش می‌رود: از غشای تحریک‌پذیر و سیستم‌های دینامیکی، به مدل‌های نورون (هاجکین–هاکسلی و مدل‌های ساده‌شده)، شبکه‌های نورون‌های اسپایکی، مدل‌های جمعیتی و توده عصبی، و سرانجام ابزارهای تحلیلی و پروژه‌های یکپارچه. هدف، فراهم‌کردن متنی روان و خودبسنده برای دانشجویان و پژوهشگران فارسی‌زبان است.

## مشارکت

این کتاب در حال توسعه است و بازخورد شما به بهتر شدن آن کمک می‌کند. برای پیشنهاد فصل یا بخش جدید، گزارش اشکال، یا هر پیشنهاد دیگری، یک [ایشو در گیت‌هاب](https://github.com/SaeedTaghavi/computational-neuroscience/issues) باز کنید یا به نشانی <saeed.taghavi.v@gmail.com> ایمیل بزنید.

## ساخت محلی

این کتاب با [MkDocs](https://www.mkdocs.org/) و تم [Material](https://squidfunk.github.io/mkdocs-material/) ساخته شده است.

```bash
pip install mkdocs-material
cd computational-neuroscience-book
mkdocs serve        # preview at http://127.0.0.1:8000
```

## انتشار

```bash
cd computational-neuroscience-book
git add -A && git commit -m "describe what changed" && git push   # save source to main
mkdocs gh-deploy --force                                          # build and publish to gh-pages
```

سایتِ منتشرشده حدود یک تا دو دقیقه پس از `gh-deploy` به‌روز می‌شود.

### نکته‌ها

- فایلِ `docs/CNAME` باید شاملِ `computational-neuroscience.ir` باشد تا دامنهٔ اختصاصی در هر انتشار حفظ شود؛ آن را حذف نکنید.
- اگر سایت پس از انتشار همچنان نسخهٔ قدیمی را نشان داد، احتمالاً کشِ Cloudflare است: از داشبوردِ Cloudflare مسیرِ Caching → Configuration → Purge Everything را اجرا کنید.
- برای `git push` از نامِ کاربریِ گیت‌هاب (`SaeedTaghavi`) و یک Personal Access Token (به‌جای رمز عبور) استفاده کنید.

## مجوز

برای جزئیاتِ مجوز با نگارنده در تماس باشید.