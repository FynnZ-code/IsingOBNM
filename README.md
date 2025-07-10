# IsingOBNM

## Simulationswerkzeug für Spin-Ketten-Modelle (Ising Modell)

Ermöglicht die Simulation und Untersuchung einer Quantenmechanischen Spin Kette im Hinblick auf ihre Quanten-Thermodynamischen Eigenschaften.
Implementiert die Operator-Basis-Norm-Methode zur Berechnung von ETH Indikatoren und Vergleicht diese mit klassischen Methoden.

## Usage

Verfügbare Methoden für -M / --Method:

MethodBasic                  | Einfache Simulation mit Magnetfeld an Position q 

ClassicalMethod              | Klassische Approximation im Energiefenster um E0

ClassicalGaussian            | Klassische Methode mit Gauß-Filter, gesteuert über Delta1/2

ClassicalFit                 | Klassische Methode mit anpassbarer Filterform (f), Tiefe (D)

Gamma                        | Simulation zur Gamma-Berechnung bei gegebener Energie

ChaosIndicator               | Berechnung eines Chaos-Indikators für die Kette

ClassicalETHDiag             | Klassische ETH-Diagonalisierung im Energiefenster

ClassicalEnergyGapHistogram  | Histogramm der Energieabstände für gegebene Länge

Beispiel:
python main.py -M ClassicalGaussian -L 6 --E0 0 --Delta 1
