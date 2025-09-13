# AlgCurve: rasterizace rovinných algebraických křivek

Program renderuje implicitní algebraické křivky v rovině a přechody mezi nimi.

## 1. Funkcionalita

- Vstup: dva polynomy, `f(x,y)` a `g(x,y)`.
- Výstup: MP4 video vizualizující přechod mezi množinami řešení `f` a `g`.

## 2. Build & run

### 2.1 Build
Spusťte `cmake --build --preset [preset]`.
Seznam možných konfigurací lze získat příkazem `cmake --list-presets`.
Pro build ve Visual Studiu se podívejte [zde](ttps://learn.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio).


#### Pokročilé build options
- `INCLUDE_DEBUG_SYMBOLS`: Zahrnout ladicí symboly (`-g3`).
- `ENABLE_SANITIZERS`: Povolit runtime checks na undefined behavior a memory safety (`-fsanitize=address,undefined`).
- `FRAME_POINTERS`: Vynutit frame pointers (`-fno-omit-frame-pointer`).

### 2.2 Run
Program přečte dva výrazy ze standardního vstupu (oddělené newline):
```bash
echo "(x^2+y^2-1) \n (x^2+y^2)^3 - 4x^2y^2" | ./algcurve
```
a vytvoří soubor `output_video.mp4` v `$(PWD)`. Dočasné snímky se ukládají do `$(PWD)/intermediate_images` a po úspěšném běhu se smažou. Pokud `ffmpeg` chybí, obrázky zůstanou a video se nevytvoří.

Knihovna křivek je v `lib/`. Příklad:
```bash
cat lib/watts_curve lib/bicorn | ./algcurve
```

#### Windows Powershell
```powershell
"(x^2+y^2-1)`r`n(x^2+y^2)^3 - 4x^2y^2" | .\algcurve
cat lib\watts_curve, lib\bicorn | .\algcurve
```

### 2.3. Syntax vstupu
Vstup je zapsán běžnou aritmetickou notací. Povolené symboly jsou proměnné (`x`, `y`), operátory `+`, `-`, `*`, `^`, závorky `(` a `)`, kde mocnina `^` vyžaduje nezáporný celý exponent, konstanty jsou zapsány v desetinném rozvoji. Mezi za sebou jdoucími konstantami a proměnnými je vloženo implicitní násobení. Příklady: `x^2 + y^2 - 1`, `(8x^3y^4 -1)^5+ 3`.

## 4. Vnitřní chod
- Na počítači nelze zpracovat všechny body intervalu [0, 1]. Program proto automaticky vybírá reprezentativní body (viz `main.cpp, image_difference`) tak, aby rozdíly mezi snímky nebyly příliš velké. Mezi dva body s velkým rozdílem (rozdíl mezi body `a, b ∈ [0, 1]` je `image_difference(render(lerp(a, f, g)), render(lerp(b, f, g))))`) se opakovaně vloží další bod, dokud všechny rozdíly neklesnou pod zadanou mez.
- Renderování křivky: Program používá test, který dokáže bezpečně vyloučit možnost existence kořenu v boxu `[-δ, δ]^2`. Pokud box může obsahovat kořen, je rozdělen na čtvrtiny a test opakován. Program se zastaví ve chvíli, kdy je box veliký jako jeden pixel. 
- Paralelizace: vlastní "threadpool" (detaily viz. `thread.hpp`) a paralelní map/for.

### 4.1 Matematické detaily vnitřního chodu
Program je založen na článku [An accurate algorithm for rasterizing algebraic curves](https://dl.acm.org/doi/pdf/10.1145/164360.164427). Výchozí test na možnou existenci kořenu poskytuje správný odhad v okolí regulárních kořenů. 
Ačkoliv je snadný na vyhodnocení, špatně odhaduje vzdálenost od kořenů s vysokou násobností (důsledkem je, že pokud má kořen násobnost `d`, vyrenderuje se jako kružnice poloměru `d`). 
Singulární body vznikají například v situaci, kde lze polynom `p` rozložit na součin `p = qr`, a existuje bod `x` t.ž. `q(x) = r(x) = 0`. Potom je bod `x` singulární kořen polynomu `p`.
V knihovně křivek lze nalézt singulární body např. v `rose`, `lines`, `circles` nebo `bifolium`.

Volitelný „desingularized“ test (v `paper.hpp` lze aktivovat odkomentováním `#define DESINGULARIZED_ROOT_CHECK`) jsem se pokusil implementovat dle článku, ale neposkytuje lepší výsledky než výchozí test.

### 4.2 Posun polynomu (a.k.a. translating the origin, recentering)
Algoritmus pro kontrolu možné existence kořene funguje pouze kolem počátku. Abychom zkontrolovali
libovolný bod `(dx, dy)` (pro polynom `p`), potřebujeme určit koeficienty posunutého polynomu `q(x, y) = p(x - dx, y - dy)`. 
To uděláme tak, že polynom `p` rozšíříme o proměnné `dx` a `dy` (nic nezměníme, ale začneme o něm uvažovat jako o prvku `ℝ[x, y, dx, dy]`), provedeme substituci `x → x - dx, y → y - dy` a následně vytkneme mocniny dx a dy tak, že dostaneme polynom s proměnnými x, y a koeficienty v `ℝ[dx, dy]`, tedy prvek `(ℝ[dx, dy])[x, y]`.
Tento výpočet je proveden jednou, a následně se pro každý bod přímočaře vyhodnotí koeficienty dosazením dx a dy.

### 4.3 Implementace polynomu

Polynom s koeficienty typu `T` v `NVARS` proměnných je reprezentován třemi vzájemně zaměnitelnými kontejnery a tenkou obálkou `Polynomial<T, NVARS>` sjednocující jejich rozhraní.

- `HashmapPolynomial<T, NVARS>`: řídká, mutable hashmapa `Monomial<NVARS> → T`. Jediný typ polynomu podporující algebraické operace. Vhodný pro abstraktní výpočty (viz. předchozí sekce), ale konstrukce a iterace přes koeficienty je pomalá.

- `SparseOssifiedPolynomial`: immutable „zmražená“ řídká reprezentace. Má fixní množinu přítomných monomů, jež je určena při konstrukci. Poskytuje `get_single_degree_slices()` (rozřezání podle celkového stupně), které se používá v desingularizovaném testu v `paper.cpp`.

- `DenseOssifiedPolynomial`: immutable „zmražená“ hustá mřížka koeficientů. Optimalizované vyhodnocení pro případ `NVARS == 2`.

- `Polynomial`: `std::variant` obálka sjednocující rozhraní. Funkce nepodporované aktuální variantou jsou vyřešeny převodem na `HashmapPolynomial`.

### 4.4 Tok dat
- String → List tokenů → AST → `HashmapPolynomial` - `input.cpp, tokenize, parse_tokens`, `main.cpp, parse_ast`
- Interpolace
- Offsetting - `render.cpp, to_offset_polynomial`
- Zmražení/ossifikace - vybráno řídké nebo husté úložiště `render.cpp, decide_ossified_polynomial_type`
- Opakované vyhodnocení a rozdělení bodů - `render.cpp, are_points_viable, subdivide_viable_points`
- Uložení bitmapy - `image.cpp, Image::write_bmp`

## 6. Hardcoded parametry
- `main.cpp`: 
  - `width = 2048` (rozlišení; snímek je `width × width`)
  - `plane_height = 4.0` (zobrazené rozpětí v ose `y`; osu `x` volí stejně širokou)
  - `max_distance = 0.1` (práh pro dělení intervalů)
  - `num_end_reps = 12` (počet duplikovaných koncových snímků kvůli pauze)
- Video: `ffmpeg -framerate 30 -c:v libx264 -pix_fmt yuv420p -crf 23`.
