# Smoothed Particle Hydrodynamics (SPH)

Dieses Projekt wurde während des Förderstipendiums “Simulierte-Welten” des Höchstleistungsrechenzentrums Stuttgart (HLRS) an der Universität Stuttgart entwickelt.

## Inhalt

### Beschreibung

Dieses Projekt ist eine Python-Implementierung der Smoothed Particle Hydrodynamics (SPH) Methode, die zur Simulation von Flüssigkeiten verwendet wird. SPH ist eine numerische Methode, die kontinuierliche Medien wie Flüssigkeiten oder Gase durch eine diskrete Anzahl von Partikeln beschreibt. Diese Partikel bewegen sich gemäß physikalischer Gesetze und Interaktionen, was es ermöglicht, komplexe Flüssigkeitsbewegungen zu simulieren, wie z.B. Strömungen, Wellen oder Spritzwasser.

### Merkmale der SPH-Methode
- Partikelbasierte Methode: Flüssigkeiten werden durch Partikel modelliert, welche Position, Geschwindigkeit und andere physikalische Größen tragen.
- Kernel-Funktionen: Diese werden verwendet, um die physikalischen Größen über die Partikel zu glätten, sodass kontinuierliche Felder wie Druck und Dichte entstehen.
- Interpartikel-Wechselwirkungen: Die Kräfte zwischen den Partikeln, wie Druck- oder Viskositätskräfte, steuern die Dynamik der Flüssigkeit.

### Anwendung der SPH-Methode
- Flüssigkeitssimulationen: Die SPH-Methode ist in der Lage, realistische Flüssigkeitsdynamiken zu simulieren, was sie ideal für Animationen und Computerspiele macht.
- Physikalisch realistische Strömungen: Mithilfe von SPH können Strömungen, Turbulenzen und Oberflächenspannungseffekte modelliert werden.

## Beispiele
### Simulation von vier verschiedenfarbigen Flüssigkeiten
In diesem Beispiel wird die SPH-Methode verwendet, um vier unterschiedliche Flüssigkeiten zu simulieren, jede mit einer eigenen Farbe. Diese Flüssigkeiten interagieren miteinander und zeigen typische Effekte wie Vermischung, Oberflächenspannung und Strömungsverhalten.

https://github.com/user-attachments/assets/82c15222-dbb4-40e7-a8e2-dce716526407

### Simulation mit Visualisierung des Drucks
Dieses Beispiel zeigt eine Flüssigkeitssimulation, bei der der Druck innerhalb der Flüssigkeit visualisiert wird. Der Druckverlauf wird als farbliches Mapping dargestellt, wobei Bereiche mit hohem Druck hervorgehoben werden. Diese Visualisierung hilft, die physikalischen Prozesse innerhalb der Flüssigkeit besser zu verstehen.

https://github.com/user-attachments/assets/bc0516e0-16d9-4bd8-bcd8-30dd27098ded

