# AlgCurve: rasterizace rovinných algebraických křivek

Program renderuje implicitně zadané algebraické křivky v rovině a přechody mezi nimi.

## 1. Funkcionalita

- Vstup: dva polynomy, `f(x,y)` a `g(x,y)`.
- Výstup: MP4 video vizualizující přechod mezi množinami řešení `f` a `g`.

## 2. Build & run

### 2.1 Build
Spusťte `cmake --build --preset [preset]`.
Seznam možných konfigurací lze získat příkazem `cmake --list-presets`.
Pro build ve Visual Studiu se podívejte [zde](ttps://learn.microsoft.com/en-us/cpp/build/cmake-projects-in-visual-studio).

Pokročilé build options:
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
Přechod mezi křivkami `f` a `g` je zařízen lineární interpolací  jejich definujících polynomů. Každý mezisnímek odpovídá jednomu bodu `t ∈ [0, 1]`, a vyobrazuje kořeny polynomu `lerp(t, f, g)`.

Na počítači nelze zpracovat všechny body intervalu `[0, 1]`. Program proto automaticky vybírá reprezentativní body (viz `main.cpp, image_difference`) tak, aby rozdíly mezi snímky nebyly příliš velké. Rozdíl mezi body `t, s ∈ [0, 1]` je `image_difference(render_image(lerp(t, f, g)), render_image(lerp(s, f, g))))`. Mezi dva body s velkým rozdílem  se opakovaně vloží další bod, dokud všechny rozdíly neklesnou pod zadanou mez.

Program používá test, který dokáže bezpečně vyloučit možnost existence kořenu v boxu velikosti `[-δ, δ]^2` kolem počátku. Pokud box může obsahovat kořen, je rozdělen na čtvrtiny (`δ[n+1]=δ[n]/2`) a test opakován. Program se zastaví ve chvíli, kdy je box veliký jako jeden pixel. Aby šlo použít test kolem počátku pro libovolný bod, je potřeba posunout polynom.

Problém je [vysoce paralelní](https://en.wikipedia.org/wiki/Embarrassingly_parallel), jelikož renderování jednotlivých obrázků ani kontrola jednotlivých bodů na sobě nijak nezávisí. Program používá [paralelní for/map](https://en.wikipedia.org/wiki/Map_(parallel_pattern)) (viz. `thread.hpp`).

### 4.1 Matematické detaily vnitřního chodu
Program je založen na článku [An accurate algorithm for rasterizing algebraic curves](https://dl.acm.org/doi/pdf/10.1145/164360.164427). Výchozí test na možnou existenci kořenu poskytuje správný odhad v okolí regulárních kořenů. 

> **Co je to singulární kořen polynomu**
> Každému kořenu `v` polynomu `p` můžeme přiřadit tzv. násobnost, značenou `m(v, p)`. Násobnost kořene je minimální index `h` t. ž. existuje `h`-tá parciální derivace, jež není v bodě `v` rovna nule.
> Ekvivalentně (pro případ dvou proměnných) to je nejvyšší index `h` t.ž. existují `a+b=h`, a `(x-v_1)^a * (y-v_2)^b` dělí `p`.

Ačkoliv je snadný na vyhodnocení, špatně odhaduje vzdálenost od singulárních kořenů (důsledkem je, že pokud má kořen násobnost `d`, vyrenderuje se jako kružnice (pixelového) poloměru `d`). 
Singulární body vznikají například v situaci, kde lze polynom `p` rozložit na součin `p = qr`, a existuje bod `x` t.ž. `q(x) = r(x) = 0`. Potom je bod `x` singulární kořen polynomu `p`.
V knihovně křivek lze nalézt singulární body např. v `rose`, `lines`, `circles` nebo `bifolium`.

Volitelný „desingularized“ test (v `paper.hpp` lze aktivovat odkomentováním `#define DESINGULARIZED_ROOT_CHECK`) jsem se pokusil implementovat dle článku, ale neposkytuje lepší výsledky než výchozí test.

### 4.2 Posun polynomu (a.k.a. translating the origin, recentering)
Algoritmus pro kontrolu možné existence kořene funguje pouze kolem počátku. Abychom zkontrolovali
libovolný bod `(dx, dy)` (pro polynom `p`), potřebujeme určit koeficienty posunutého polynomu `q(x, y) = p(x - dx, y - dy)`. 
> **Okruh polynomů R[x_1, x_2, ..., x_n]**
> Pro libovolný okruh **R** a celé číslo **n** definujeme okruh 
> polynomů nad **R** v **n** proměnných, značeno **R[x_1, ..., x_n]**, jako množinu konečných posloupností **Σc[i_1, ..., i_n] * x_1^i_1 * x_2^i_2 * ... * x_n^i_n**.

To uděláme tak, že polynom `p` rozšíříme o proměnné `dx` a `dy` (nic nezměníme, ale začneme o něm uvažovat jako o prvku `ℝ[x, y, dx, dy]`), provedeme substituci `x → x - dx, y → y - dy` a následně vytkneme mocniny dx a dy tak, že dostaneme polynom s proměnnými x, y a koeficienty v `ℝ[dx, dy]`, tedy prvek `(ℝ[dx, dy])[x, y]`.
> Příklad substituce: z polynomu **3x^2^.dx+4x.y.dy+dx+2x.y^2^+7+9x^2^.dx.dy^5^+5x^2^.dx.dy**
> uděláme polynom **(3dx+9dy^5^+5dx.dy)x^2^+(4dy)x.y+(2)x.y^2^+(dx+7)1**.

Tento výpočet je proveden jednou, a následně se pro každý bod přímočaře vyhodnotí koeficienty dosazením dx a dy.

### 4.3 Implementace polynomu

Polynom s koeficienty typu `T` v `NVARS` proměnných je reprezentován třemi vzájemně zaměnitelnými kontejnery a tenkou obálkou `Polynomial<T, NVARS>` sjednocující jejich rozhraní.

- `HashmapPolynomial<T, NVARS>`: řídká mutable hashmapa `Monomial<NVARS> → T`. Jediný typ polynomu podporující algebraické operace. Vhodný pro abstraktní výpočty (viz. předchozí sekce), ale konstrukce a iterace přes koeficienty je pomalá.

- `SparseOssifiedPolynomial`: řídká immutable reprezentace. Má fixní množinu přítomných monomů, jež je určena při konstrukci. Poskytuje `get_single_degree_slices()` (rozřezání podle celkového stupně), které se používá v desingularizovaném testu v `paper.cpp`.

- `DenseOssifiedPolynomial`: hustá immutable reprezentace. Má fixní stupně jednotlivých proměnných. Optimalizované vyhodnocení pro případ `NVARS == 2`.

- `Polynomial`: `std::variant` obálka sjednocující rozhraní. Funkce nepodporované aktuální variantou jsou vyřešeny převodem na `HashmapPolynomial`.

### 4.4 Tok dat
- String → List tokenů → AST → `HashmapPolynomial` - `input.cpp, tokenize, parse_tokens`, `main.cpp, parse_ast`
- Interpolace
- Offsetting - `render.cpp, to_offset_polynomial`
- Zmražení/ossifikace - vybráno řídké nebo husté úložiště `render.cpp, decide_ossified_polynomial_type`
- Opakované vyhodnocení a rozdělení bodů - `render.cpp, are_points_viable, subdivide_viable_points`
- Uložení bitmapy - `image.cpp, Image::write_bmp`

## 6. Hardcoded parametry
Některé parametry programu jsou zafixované při kompilaci.
- `width = 2048` (rozlišení; snímek je `width × width`)
- `plane_height = 4.0` (zobrazené rozpětí v ose `y`; osu `x` volí stejně širokou)
- `max_distance = 0.1` (práh pro dělení intervalů)
- `num_end_reps = 12` (počet duplikovaných koncových snímků kvůli pauze)
- Video: `ffmpeg -framerate 30 -c:v libx264 -pix_fmt yuv420p -crf 23`.
