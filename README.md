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

### 4.3 Implementační detaily
```
┌─────────────┐   ┌──────────────┐                            
│             │   │              ├──────────────────────┐     
│  heaparray  ├──►│  polynomial  │                      │     
│             │   │              ├────────┐             │     
└─────────────┘   └───────┬──────┘        │             │     
                          │               │             │     
                          ▼               ▼             ▼     
┌─────────┐       ┌──────────┐       ┌─────────┐   ┌─────────┐
│         │       │          │       │         │   │         │
│  image  ├──────►│  render  │◄──────┤  paper  │   │  input  │
│         │       │          │       │         │   │         │
└─────────┘       └──┬───────┘       └─────────┘   └────┬────┘
                     │    ▲                             │     
                     │    │          ┌──────────┐       │     
┌────────┐◄──────────┘    └──────────┤          │       │     
│        │                           │  thread  │       │     
│  main  │◄──────────────────────────┤          │       │     
│        │                           └──────────┘       │     
└────────┘◄─────────────────────────────────────────────┘     
```
#### 4.3.1 Zarovnané úložiště, `core.hpp, SimdHeapArray`
Třída `SimdHeapArray` reprezentuje pole fixní velikosti konzervativně zarovnané pro použítí instukcí SIMD. Zarovnané alokace fungují jinak na Windows, což vyžaduje čarování s makry.
#### 4.3.2 Reprezentace polynomu, `polynomial.hpp`

Polynom s koeficienty typu `T` v `NVARS` proměnných je reprezentován třemi vzájemně zaměnitelnými kontejnery a tenkou obálkou `Polynomial<T, NVARS>` sjednocující jejich rozhraní.

- `HashmapPolynomial<T, NVARS>`: řídká mutable hashmapa `Monomial<NVARS> → T`. Jediný typ polynomu podporující algebraické operace. Vhodný pro abstraktní výpočty (viz. předchozí sekce), ale konstrukce a iterace přes koeficienty je pomalá.

- `SparseOssifiedPolynomial`: řídká immutable reprezentace. Má fixní množinu přítomných monomů, jež je určena při konstrukci. Poskytuje `get_single_degree_slices()` (rozřezání podle celkového stupně), které se používá v desingularizovaném testu v `paper.cpp`.

- `DenseOssifiedPolynomial`: hustá immutable reprezentace. Má fixní stupně jednotlivých proměnných. Optimalizované vyhodnocení pro případ `NVARS == 2`.

- `Polynomial`: `std::variant` obálka sjednocující rozhraní. Funkce nepodporované aktuální variantou jsou vyřešeny převodem na `HashmapPolynomial`.

#### 4.3.3 Reprezentace obrázku, `image.hpp`
Program produkuje bitmapy s přednastavenou paletou, jeden pixel je reprezentován jako jeden byte. Třída `Image` používá hustou reprezentaci obrázku, ačkoliv pro čtvercový obrázek velikosti `n` jednorozměrná křivka zabere pouze `O(n)` pixelů.

#### 4.3.4 Multithreading, `thread.hpp`
Každé vlákno má přiřazené, kolik pracujících vláken smí spustit. Když se zavolá `parallel_for` s méně jednotkami práce, než je dostupných vláken, "přebytečná" vlákna se rozdělí nově vzniklým vláknům. To se typicky děje při renderování malého množství obrázků v jendom batchi, kde poté dojde k paralelizaci kontroly kořenů.

#### 4.3.5 Vstup, `input.hpp`
Jednoduchý parser vstupu. Provede tokenizaci, převede seznam tokenů na AST a ten následně vyhodnotí na polynom. Tokeny ukládají, kde se ve vstupním textu nacházejí, což umožňuje generovat základní chybové zprávy.

#### 4.3.6 Algoritmus, `paper.hpp`
Zde jsou implementovány algoritmy pro možnou existenci kořene kolem počátku ze článku An accurate algorithm for rasterizing algebraic curves. Proces offsettingu je odložen do `render.cpp`.

#### 4.3.7 Rasterizace, `render.hpp`
Rasterizace je řízená dvojicí funkcí `render_image` a `render_images`. Pro rychlé vyhodnocování posunutého polynomu `p(x-dx, y-dy)` se připraví tabulky mocnin souřadnic mřížky do struktury `PreparedLattices`. Každá granularita má vlastní pravidelnou mřížku bodů se středem v každé buňce. Šířka obrázku musí být mocninou dvou, aby mřížka přesně seděla na pixely.

`render_image` implementuje dělení a ověřování kandidátních buněk: z počátečního bodu na granularitě 0 se opakovaně zavolá test možného kořene (`may_have_root`) nad boxem kolem bodu, a „průchozí“ body se rozdělí na čtyři podbody (`subdivide_viable_points`). "Poloměr" boxu `δ` je polovina šířky buňky. Smyčka končí při dosažení granularitě odpovídající rozlišení (`log2(width)`), kdy každý zbývající bod odpovídá jednomu pixelu, jenž se obarví na `WHITE`.

`render_images` vyrenderuje sadu snímků pro zadané hodnoty `λ`. Každý snímek je reprezentován polynomem `p = lerp(λ, p1, p2)`. Ten se nejprve převede na posunutelnou reprezentaci `OffsetPolynomial` přes `to_offset_polynomial`, a následně se zavolá `render_image` s předpočítanými mřížkami.

#### 4.3.8 Jádro, `main.cpp`
`main` koordinuje celý běh. Načte dva výrazy ze standardního vstupu, parsne je na `HashmapPolynomial<double,2>` a z jejich stupně odvodí parametry pro `PreparedLattices`.

Pro adaptivní volbu klíčových snímků se nejdřív vyrenderují okrajové hodnoty `t=0` a `t=1` a uloží do mapy `interpolation_steps: t → Image`. Následně program udržuje frontu intervalů `[l, r]`, pro každý spočítá metrikou `image_difference(l, r)` rozdíl obrázků (paralelně), a pokud je nad prahem `max_distance`, vloží střed `mid` a rozdělí interval. Všechny nově vzniklé středy `mid` se vyrenderují v batchích přes `render_images` a doplní do mapy. Smyčka končí, když žádný interval nepřekračuje práh.

Po konvergenci se do každého snímku doplní progress bar (podle parametru `t`), připraví se dopředná i zpětná sekvence s krátkou pauzou na koncích (`num_end_reps`) a snímky se uloží jako `intermediate_images/zzz%04d.bmp` (paralelně). Pokud je v systému dostupné `ffmpeg`, vytvoří se z nich `output_video.mp4`; v opačném případě program ponechá BMP soubory a uživatele informuje. Dočasné soubory se po úspěchu uklidí. 

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

## 7. Rekapitulace
Na myšlenku renderování algebraických křivek jsem přišel při čtení úvodní sekce knihy [A guide to plane algebraic curves](https://agorism.dev/book/math/curve/guide-to-plane-algebraic-curves_keith-kendig.pdf). Zde využitý článek nabídl snadné řešení problému rasterizace, ale algoritmus pro práci se singulárními body se ukázal být náročnější, než jsem čekal, a nezvládl jsem ho správně implementovat. 
Myšlenka videa přechodu také pochází z knihy o křivkách, a poskytla zajímavý problém ukládání snímků obrázků na renderování - program dokáže vygenerovat obří množství dat (což se dá vyřešit řídkým kódováním obrázků), které je i pomalé uložit na disk, což mě nečekaně zavedlo k RLE kódování bitmap.
Taky jsem si znovu připomněl krásy C++ konstruktorů a memory managementu, a neplánovaný vývoj softwaru na 3 různých platformách mě dovedl k prvnímu použití CMake. Nebylo to tak hrozné, jak jsem slyšel.
Celkově jsem se na práci s C++ moc netěšil, ale i paralelizace se ukázala být s použitím jednoduchého paralelního foru zvládnutelná. ASan, UBSan a debugger byly nedocenitelné nástroje.