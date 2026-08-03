_LOGO = """
╔═════════════════════════════════════════════════════════════════════════╗
║                                                                         ║
║    ██╗███████╗ ██████╗  ██████╗ ██╗   ██╗ █████╗ ███╗   ██╗████████╗    ║
║    ██║██╔════╝██╔═══██╗██╔═══██╗██║   ██║██╔══██╗████╗  ██║╚══██╔══╝    ║
║    ██║███████╗██║   ██║██║   ██║██║   ██║███████║██╔██╗ ██║   ██║       ║
║    ██║╚════██║██║   ██║██║▄▄ ██║██║   ██║██╔══██║██║╚██╗██║   ██║       ║
║    ██║███████║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║██║ ╚████║   ██║       ║
║    ╚═╝╚══════╝ ╚═════╝  ╚══▀▀═╝  ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═══╝   ╚═╝       ║
║                                                                         ║
║    BaseCode IsoQuant Pipeline · Powered by IsoQuant                     ║
{ver_line}
{ruo_block}
║    © 2026 Basic Genomics AB · All rights reserved                       ║
╚═════════════════════════════════════════════════════════════════════════╝"""

_RUO = ["Basic Genomics' software, products and services are",
        "for research use only and not for use in diagnostic procedures."]

def _box(text=""):
    return "║" + (("    " + text) if text else "").ljust(73) + "║"

def build_banner(version, codename=""):
    ver_line = _box(" ".join(filter(None, ["v" + version, codename, "Release"])))
    blank = _box()
    ruo_block = "\n".join([blank] + [_box(t) for t in _RUO])
    return _LOGO.format(ver_line=ver_line, ruo_block=ruo_block)
