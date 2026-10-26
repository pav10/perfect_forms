"""The shared cypari2 instance (one per process)."""
import cypari2

pari = cypari2.Pari()
try:
    pari.allocatemem(2 ** 28, 2 ** 32, silent=True)
except TypeError:  # older cypari2
    pari.allocatemem(2 ** 28, 2 ** 32)
except Exception:
    pass
pari.setrand(1)
