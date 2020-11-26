\LoadClass[10pt]{article}

\usepackage{graphicx}

\ProvidesClass{utiarr}[2004/12/01 v1.01]


\let\@email\@empty
\newcommand{\email}[1]{\gdef\@email{#1}}


\let\@rrno\@empty
\newcommand{\rrno}[1]{\gdef\@rrno{#1}}


\renewcommand{\maketitle}
{
\begin{titlepage}                       % "titlepage" - textheight 620pt
\vspace*{-40pt}                         % eliminate top margin setting
\includegraphics[width=40mm]{logo} 
\hspace*{4.2cm}
\begin{minipage}{8.5cm}
\sf
Akademie v\v{e}d \v{C}esk\'e republiky\\
\'{U}stav teorie informace a automatizace\\[0.4cm]
Academy of Sciences of the Czech Republic\\
Institute of Information Theory and Automation
\end{minipage}
\vspace*{1cm}
\vfill
\begin{center}
{\LARGE\sf RESEARCH REPORT}\\[0.5cm]
\vfill
\fbox{\parbox{13cm}{
\vspace*{1.3cm}
\vfill
\begin{center}
\large\sc
%%%%%%%% zde napsat jm\'eno autora(\accent23u) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\@author \\

\vspace{0.6cm} 

\Large\bf \addtolength{\baselineskip}{0.1cm}        % 

\@title 

\end{center}
\vfill
\vspace*{1.3cm}
}}                 % end of fbox, parbox
%%%%%%%% n\'asleduje pr\'azdn\'y \v{r}\'adek - ukon\v{c}en\'{\i} odstavce pro dal\v{s}\'{\i} "vspace"

\vspace*{1cm}
\fbox{\parbox{13cm}{
\vspace*{0.5cm}
\hspace*{1cm}\Large\sf No.
%%%%%%%% zde napsat \v{c}\'{\i}slo zpr\'avy, nap\v{r}. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\@rrno
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\hfill
%
%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%
\@date
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\hspace*{0.7cm}
\vspace*{0.5cm}
}}   % end of fbox, parbox
\vfill
\vspace*{1.5cm}
{\sf \'{U}TIA  AV \v{C}R, P. O. Box 18, 182 08 Prague, Czech Republic\\
     Telex: 122018 atom c,~~Fax: (+420) 266 053 111\\
     E-mail: utia@utia.cas.cz}
\end{center}
\end{titlepage}
%\setlength{\evensidemargin}{38.5pt}    % active from this page - optional
%\setlength{\oddsidemargin}{38.5pt}     % active from this page - optional


\newpage
\thispagestyle{empty}

% vnit\v{r}n\'{\i} strana p\v{r}edn\'{\i}ch desek

This report constitutes an unrefereed manuscript which is intended
to be submitted for publication. Any opinions and conclusions
expressed in this report are those of the author(s) and do not
necessarily represent the views of the Institute.


\newpage
}




\topmargin 40pt               % default value 27pt 

\evensidemargin 18pt          % default value 39.5pt 

\oddsidemargin  18pt          % default value 39.5pt 

\parskip 3pt plus 1pt         % default value 0pt plus 1pt 

\headheight 0pt               % default value 12pt 

\headsep 0pt                  % default value 25pt 

\footskip  40pt               % default value 30pt 

\textheight 620pt             % default value 536.5pt 

\textwidth 430pt              % default value 390pt
     % these values for "titlepage"; can be changed to other values later
%\def\baselinestretch{1.0}    % default value {1} - single-spaced

\def\insertplot#1#2#3{\par  % graphs in default (TEXINPUT) or explicitly
  \hbox{                                       % configurated directory
    \hskip #3               % hskip from left margin to picture left rand
    \vbox to #2 {           % vertical space dedicated for picture
      { \special{em:graph #1}}
       \vfil
      }  }  }                % insertplot
%\pagestyle{empty}
