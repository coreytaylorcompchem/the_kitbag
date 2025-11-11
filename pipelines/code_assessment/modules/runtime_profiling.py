import cProfile
import pstats
import io
import subprocess

def profile_pipeline(entrypoint_script):
    profile = cProfile.Profile()
    profile.enable()

    subprocess.run(["python", entrypoint_script])

    profile.disable()

    s = io.StringIO()
    ps = pstats.Stats(profile, stream=s).sort_stats("cumulative")
    ps.print_stats(20)  # top 20 hotspots

    return {"profile_summary": s.getvalue()}
