\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{bumazhka}[2013/12/13 Burecraucy document]
\def\my@fontsize{12pt}\def\my@baseclass{article}
\DeclareOption{10pt}{\def\my@fontsize{14pt}}
\DeclareOption{14pt}{\def\my@fontsize{14pt}\def\my@baseclass{extarticle}}
\ProcessOptions\relax
\LoadClass[a4paper,\my@fontsize]{\my@baseclass}

\usepackage[a4paper, top=1cm, bottom=1cm, left=2cm, right=2cm]{geometry}
\usepackage{amsmath}
\usepackage{dashrule}
\usepackage{fontspec}
\usepackage{polyglossia}
\usepackage{hyperref}
\usepackage{calc}

\usepackage{array}
\usepackage{xkeyval}

\makeatletter
\define@key{lines_env}{w}{\def\mywidth{#1}}
\define@choicekey*{lines_env}{align}[\val\nr]{left,center,right}{%
  \ifcase\nr\relax
    \def\myalignment{flushleft}
  \or
    \def\myalignment{center}
  \or
    \def\myalignment{flushright}
  \fi
}
\makeatother

\newenvironment{lines}[1][]
{\setkeys{lines_env}{#1}
\newcolumntype{E}{>{}p{\mywidth}<{\hrulefill}}
\begin{\myalignment}
\begin{tabular}[h]{E}} 
{\end{tabular}\end{\myalignment}}

\usepackage{xstring}
\usepackage{collcell}
\newcommand{\MyCellRenderer}[1]{#1}%
\newcolumntype{F}[1]{>{\collectcell\MyCellRenderer}p{#1}<{\endcollectcell}}%
\def\myspec{}
\def\doword#1{\expandafter\def\expandafter\myspec\expandafter{\myspec F{#1\textwidth}}}
\def\endtestwords{}
\def\readcolspec#1{%
    \begingroup
    \edef\@tempa{#1;}%
    \expandafter\endgroup%
    \expandafter\readwords\@tempa\relax
}
\def\readwords#1;#2\relax{%
      \doword{#1}%  #1 = substr, #2 = rest of string
      \begingroup
      \ifx\relax#2\relax  % is #2 empty?
         \def\next{\endgroup\endtestwords}% your own end-macro if required
      \else
         \def\next{\endgroup\readwords#2\relax}%
      \fi
      \next
}

\makeatletter
\define@key{cells_env}{w}{\def\MyCellsWidth{#1}}
\define@key{cells_env}{spec}{\readcolspec{#1}}
\makeatother

\newenvironment{cells}[1][]
{\setkeys{cells_env}{#1}\renewcommand*{\arraystretch}{2.0}
\expandafter\tabular\expandafter{\myspec}
}{
\endtabular
}


\setdefaultlanguage{russian}
\setmainfont[Mapping=tex-text]{CMU Serif}
\newcommand\chkbox[1]{\fbox{\begin{minipage}[c][0.3cm]{0.3cm}\quad\end{minipage}}\:\parbox[t]{\linewidth - 0.3cm}{#1}
}
\newcommand\chkboxx[1]{\begin{minipage}[c][0.3cm]{0.3cm}\quad\end{minipage}\:\parbox[t]{\linewidth - 0.3cm}{#1}
}

\setlength{\tabcolsep}{2pt}
\setlength{\parindent}{0pt}
