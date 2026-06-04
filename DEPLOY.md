# Deploy notes (private)

Personal reference for publishing the site. Not meant for readers.

## Publish

```bash
cd computational-neuroscience-book
git add -A && git commit -m "describe changes" && git push
mkdocs gh-deploy --force
```

Live site updates ~1–2 min after `gh-deploy`.

## Gotchas

- `docs/CNAME` must contain the custom domain so it survives every deploy. Don't delete it.
- If the live site looks stale after deploying, purge the CDN cache (Cloudflare → Caching → Purge Everything).
- Push auth uses a Personal Access Token, not a password.

## Local preview

```bash
cd computational-neuroscience-book
mkdocs serve        # http://127.0.0.1:8000
```
