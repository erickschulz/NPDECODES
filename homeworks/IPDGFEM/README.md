## Homework IPDGFEM for NumPDE course

* Bearbeite nur die files in `developers/IPDGFEM/`.
* Neue files und libraries ins dependencies.cmake schreiben!
* Um die templates von mastersolution zu unterscheiden, schreibe #if SOLUTION tags.
* Als Beispiel nimm BurgersEquation/
* Der erste build kann es sehr lange dauern, da viele libraries installiert werden müssen!

Um deine executables zu generieren, navigiere im Terminal ins höchste Vezeichnis des Repositories und dann:
```
mkdir build
cd build
cmake -DHOMEWORKS=OFF -DLECTURECODES=OFF -DMYSOLUTION=OFF ..
cd developers/IPDGFEM/
make
```
Dies generiert die executables `IPDGFEM_mastersolution` und `IPDGFEM_test_mastersolution` im aktuellen Verzeichnis, d.h. `build/developers/IPDGFEM/`.
