# Dependencies required (all present in pyproject.toml):
#   fastapi>=0.100, uvicorn[standard]>=0.20, python-multipart>=0.0.5
"""PairMap2 Web UI — FastAPI application entry point."""
from __future__ import annotations

import logging
import os
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

logger = logging.getLogger(__name__)

from . import executor, job_store, map_store
from .routes import health, jobs, pair, map as map_routes

app = FastAPI(title="PairMap")


@app.exception_handler(Exception)
async def unhandled_exception_handler(request: Request, exc: Exception):
    logger.exception("Unhandled exception for %s %s", request.method, request.url)
    return JSONResponse(status_code=500, content={"detail": str(exc) or "An unexpected error occurred."})

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://pairmap.yumizsui.com",
        "http://localhost:8000",
        "http://localhost:3000",
    ],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Health + PairMap2 jobs (existing)
app.include_router(health.router)
app.include_router(jobs.router, prefix="/api")

# New: Pair Mode (sync) and Map Mode v1
app.include_router(pair.router, prefix="/api")
app.include_router(map_routes.router, prefix="/api/map")


@app.on_event("startup")
async def startup() -> None:
    job_store.init_db()
    map_store.init_db()
    purged = map_store.purge_old_jobs()
    if purged:
        import logging
        logging.getLogger(__name__).info("Purged %d old map jobs", len(purged))
    executor.init()


@app.on_event("shutdown")
async def shutdown() -> None:
    executor.shutdown()


# Serve frontend SPA
frontend_path = Path(__file__).parent.parent / "frontend"
if frontend_path.exists():
    app.mount(
        "/",
        StaticFiles(directory=str(frontend_path), html=True),
        name="frontend",
    )
