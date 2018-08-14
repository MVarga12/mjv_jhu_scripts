## Instructions for using Minion Pro in LaTeX
  - Minion Pro needs to be mined from Adobe Acrobat Reader (or any other Adobe software)
  - Can be included in LaTeX with `\usepackage{MinionPro}`
    - I use `\usepackage[smallfamily,noopticals,lf,italicgreek,loosequotes]{MinionPro}`

### Minion Pro for Math
  - Need to explicitly set Minion Pro to be the math font:
    ```
    \usepackage{unicode-math}
    \setmainfont[Numbers = OldStyle,Ligatures = TeX,SmallCapsFeatures = {Renderer=Basic}]{Minion Pro}
    \setmathfont{MnSymbol}
    \setmathfont[range=\mathup/{num,latin,Latin,greek,Greek}]{Minion Pro}
    \setmathfont[range=\mathbfup/{num,latin,Latin,greek,Greek}]{MinionPro-Bold}
    \setmathfont[range=\mathit/{num,latin,Latin,greek,Greek}]{MinionPro-It}
    \setmathfont[range=\mathbfit/{num,latin,Latin,greek,Greek}]{MinionPro-BoldIt}
    \setmathrm{Minion Pro}
    ```
