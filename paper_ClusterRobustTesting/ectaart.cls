%%
%% LaTeX 2e class file for the processing of LaTeX2e files
%% for the Econometric Society journal "Econometrica"
%%
%% Macros written by Vytas Statulevicius, VTeX, Lithuania
%% for the Econometric Society, U.S.A.
%% Please submit bugs or your comments to vytas@vtex.lt
%%
%% The original distribution is located at:
%% http://www.e-publications.org/ecta/support
%%
%% This class file loads standart "article.cls" with appropriate 
%% settings and then redefines layout according to "Econometrica" style
%%
%% You are free to use this style class as you see fit, provided 
%% that you do not make changes to the file. 
%% If you DO make changes, you are required to rename this file.
%%
%% It may be distributed under the terms of the LaTeX Project Public
%% License, as described in lppl.txt in the base LaTeX distribution.
%% Either version 1.0 or, at your option, any later version.
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
%%
%%
%% Bug fixes and changes:

% TeX programming: Vytas Statulevicius, VTeX, Lithuania, vytas@vtex.lt
% Requires Latex2e, ver.2000.06

\def\ectafmt@name{ectaart.cls}
\def\ectafmt@version{2006/04/11}

\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{ectaart}
   [\ectafmt@version Econometrica article class (VS)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% General options:

% Print id line at bottom of the page:
\DeclareOption{noinfoline}{\AtBeginDocument{\let\info@line\@empty}}
\DeclareOption{infoline}  {\AtBeginDocument{\let\info@line\infoline@text}}

% Put lines numbers in margins
\newif\ifnumberlines@ \numberlines@false
\DeclareOption{linenumbers}{\numberlines@true}
\DeclareOption{nolinenumbers}{\numberlines@false}

% Spacing
\DeclareOption{doublespacing}{\AtBeginDocument{\renewcommand{\baselinestretch}{1.4}\large\normalsize}}
\DeclareOption{singlespacing}{\AtBeginDocument{\renewcommand{\baselinestretch}{1.0}\large\normalsize}}

% Size - if option 10pt is used, true journal dimensions are kept
\newif\if@xiisize \@xiisizefalse

% Change in layout: put addresses as footnotes
\newif\if@addressasfootnote
\@addressasfootnotetrue
\DeclareOption{addressasfootnote}{\@addressasfootnotetrue}                  
\DeclareOption{addressatend}     {\@addressasfootnotefalse}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Options for users: draft, final

% "draft" option must be used for submission:
\newif\if@ectadraft \@ectadraftfalse
                                             
\DeclareOption{draft}{%
  \AtBeginDocument{\let\appendix@size\normalsize}
  \AtBeginDocument{\renewcommand{\baselinestretch}{1.4}\large\normalsize}
  \@xiisizetrue
  \numberlines@true
  \@ectadrafttrue
  \AtBeginDocument{\let\info@line\infoline@text}
  \PassOptionsToClass{12pt,titlepage}{article}
  \PassOptionsToPackage{draft=false}{hyperref}
  \PassOptionsToPackage{final}{graphics}
}

% "final" option must be used for posting of accepted version:
\DeclareOption{final}{%
  \@xiisizefalse
  \numberlines@false
  \PassOptionsToClass{10pt}{article}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading standart "article" class


% Set draft to default
%\ExecuteOptions{draft}

\PassOptionsToClass{fleqn,leqno,twoside}{article}
\PassOptionsToPackage{fleqn,leqno}{amsmath}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*

\LoadClass{article}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Modifications and "add-on" for article.cls starts:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initiate some info:
\def\journal@name{Submitted to Econometrica}
\def\journal@url{http://www.econometricsociety.org/}
\def\journal@issn{0012-9682}
\def\journal@id{-econ}
\def\paper@url{}
\def\info@line{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Page dimensions

\if@xiisize
  \setlength\textwidth{32pc}
  \setlength\textheight{48pc}
\else
  \setlength\textwidth{29pc}
  \setlength\textheight{43pc}
\fi
  \advance\textheight by\topskip

% make side margins equal:
  \setlength\@tempdima        {\paperwidth}
  \addtolength\@tempdima      {-\textwidth}
  \setlength\oddsidemargin    {.5\@tempdima}
  \addtolength\oddsidemargin  {-1in}
  \setlength\evensidemargin   {\oddsidemargin}
  \@settopoint\oddsidemargin
  \@settopoint\evensidemargin

% topmargin
  \setlength\topmargin{\paperheight}
  \addtolength\topmargin{-1in}
  \addtolength\topmargin{-\headheight}
  \addtolength\topmargin{-\headsep}
  \addtolength\topmargin{-\textheight}
  \addtolength\topmargin{-\footskip}     % this might be wrong!
  \addtolength\topmargin{-.5\topmargin}
  \@settopoint\topmargin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dimensions

\if@xiisize
  \def\true@parindent  {12\p@}
  \setlength\parindent {12\p@}
  \setlength\headsep   {15\p@}
  \setlength\footskip  {26\p@}

  \setlength\smallskipamount{7\p@ \@plus 1\p@ \@minus 1\p@}
  \setlength\medskipamount  {14\p@ \@plus 3\p@ \@minus 3\p@}
  \setlength\bigskipamount  {21\p@ \@plus 3\p@ \@minus 3\p@}
\else
  \def\true@parindent  {10\p@}
  \setlength\parindent {10\p@}
  \setlength\headsep   {14\p@}
  \setlength\footskip  {22\p@}

  \setlength\smallskipamount{6\p@ \@plus 1\p@ \@minus 1\p@}
  \setlength\medskipamount  {12\p@ \@plus 3\p@ \@minus 3\p@}
  \setlength\bigskipamount  {18\p@ \@plus 3\p@ \@minus 3\p@}
\fi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setattribute, getattribute, do@option@list

\def\setattribute{\@ifnextchar[\@setattribute{\@setattribute[]}}
\def\@setattribute[#1]#2#3#4{\expandafter\gdef\csname #2@#3\endcsname{#4}}
\def\getattribute#1#2{\csname #1@#2\endcsname}
\def\sep@key@value#1=#2/?/#3{\setattribute{#3}{#1}{#2}}
\def\do@option@list#1#2{%
  \@for\curr@option:={#2}\do{%
    \expandafter\sep@key@value\curr@option/?/{#1}\relax
  }%
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% newpseudoenvironment
% same as \newenvironment, but new environment do not have additional groups \bgroup \egroup
% (i.e. all definitions are not local

\let\org@begin\begin
\let\org@end\end
\def\begin#1{%
  \@ifundefined{pseudo@#1}%
    {\org@begin{#1}}{\csname pseudo@#1\endcsname[0]\relax}%
  }
\def\end#1{%
  \@ifundefined{pseudo@#1}%
    {\org@end{#1}}{\csname pseudo@#1\endcsname[1]\relax}%
  }
\def\newpseudoenvironment#1#2#3{%
  \expandafter\gdef\csname pseudo@#1\endcsname[##1]{%
     \relax\ifcase##1\relax\def\@@next@@{#2}\or\def\@@next@@{#3}\else\let\@@next@@\relax\fi\@@next@@}%
  }


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pubyear, volume, paperno

\def\pubmonth#1{\gdef\@pubmonth{#1}}%
\gdef\@pubmonth{}

\def\pubyear#1{\gdef\@pubyear{(#1) }\gdef\@copyrightyear{#1 }}
  \def\@pubyear{}
  \def\@copyrightyear{}

\def\volume#1{\gdef\@volume{Vol. #1 }}
  \gdef\@volume{}

\def\issue#1{\gdef\@issue{No. #1}}
  \gdef\@issue{}

\def\paperno#1{\gdef\@paperno{Paper no. #1 }}
  \gdef\@paperno{00} 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% firstpage, lastpage, pagerange

\newcounter{firstpage}
\newcounter{lastpage}

\def\firstpage#1{\def\@tempa{#1}\ifx\@tempa\@empty\else
  \global\c@firstpage=#1
  \global\c@lastpage=#1
  \global\c@page=#1 \ignorespaces\fi}

\def\lastpage#1{\def\@tempa{#1}\ifx\@tempa\@empty\else
  \global\c@lastpage=#1
  \ignorespaces\fi}

\def\pagerange@sep{--}

\def\set@pagerange{%
  \ifnum\c@firstpage=0%
  \else%
     \ifnum\c@firstpage=\c@lastpage%
        \gdef\@pagerange{\thefirstpage}%
     \else%
        \gdef\@pagerange{\thefirstpage\pagerange@sep\thelastpage}%
     \fi%
   \fi}

\def\@pagerange{}

\def\pagenumbering#1{%
    \gdef\thefirstpage{\csname @#1\endcsname\c@firstpage}%
    \gdef\thelastpage{\csname @#1\endcsname\c@lastpage}%
    \gdef\thepage{\csname @#1\endcsname\c@page}%
}

% hyperref redefines \pagenumbering, so we must override hyperref definition:
\let\ims@pagenumbering\pagenumbering


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% startlocaldefs, endlocaldefs

\def\startlocaldefs{\makeatletter}
\def\endlocaldefs{\makeatother}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% thanksref, thanksmark, thankslabel, thankstext
% to be safe with hyperref we will use original LaTeX definitions:
%

\def\saferef#1{\expandafter\safe@setref\csname r@#1\endcsname\@firstoftwo{#1}}
\let\safe@setref\@setref

\def\safelabel#1{\@bsphack
  \protected@write\@auxout{}%
         {\string\thanksnewlabel{#1}{{\@currentlabel}{\thepage}}}%
  \@esphack}

\let\thanksnewlabel\newlabel

% we want to use various counters:
\def\usethankscounter#1{%
  \@ifundefined{current@thankscounter}{\gdef\previous@thankscounter{#1}}{\xdef\previous@thankscounter{\current@thankscounter}}%
  \def\current@thankscounter{#1}}

\def\restorethankscounter{\xdef\current@thankscounter{\previous@thankscounter}}

\newcounter{thanks}
%\def\thethanks{\@fnsymbol\c@thanks}
\def\thethanks{\@arabic\c@thanks}
\usethankscounter{thanks}

\def\thanksmark@fmt#1{\hbox{$^{#1}$}}
\def\thanksref@sep{,}

% hooks for the hyperref:
\def\thankref@hyperlink#1{\saferef{#1thanks}}
\def\thanks@hypertarget#1{}

% Isvedame zymes
\def\thanksref{\@ifnextchar[{\@tempswatrue\@thanksref}{\@tempswafalse\@thanksref[]}}

\def\@thanksref[#1]#2{%
  \if@tempswa% []
    \thanksmark@fmt{#1}%
  \else%
    \let\@tempa\@empty%
    \thanksmark@fmt{\@for\@tempb:=#2\do{%
       \@tempa\let\@tempa\thanksref@sep%
       \edef\@tempb{\expandafter\@firstofone\@tempb\@empty}%
       \thankref@hyperlink{\@tempb}}}%
   \fi}

% Suformuojame ir isvedame zyme
\def\thanksmark{\@ifnextchar[{\@tempswatrue\@thanksmark}{\@tempswafalse\@thanksmark[]}}

\def\@thanksmark[#1]#2{%  
   \@thankslabel[#1]{#2}%
   \safelabel{#2thanks}%
   \thanksmark@fmt{\expandafter\saferef{#2thanks}\thanks@hypertarget{#2}}}

% Suformuojame tik zyme
\def\thankslabel{\@ifnextchar[{\@tempswatrue\@thankslabel}{\@tempswafalse\@thankslabel[]}}

\def\@thankslabel[#1]#2{%
  \if@tempswa% []
     \protected@edef\@currentlabel{#1}%
   \else% 
      \refstepcounter{\current@thankscounter}%
   \fi%
   \safelabel{#2thanks}}%

% Suformuojame zyme ir idedame teksta i \@thanks:
\def\thankstext{\@ifnextchar[{\@tempswatrue\@thankstext}{\@tempswafalse\@thankstext[]}}

\def\@thankstext[#1]#2#3{%
  \@thankslabel[#1]{#2}%
  \protected@xdef\@thanks{\@thanks\protect\thanks@thefnmark{#2thanks}%
  \protect\@footnotetext{\thanks@hypertarget{#2}#3}}}%

\def\thanks@thefnmark#1{\begingroup\unrestored@protected@xdef\@thefnmark{\saferef{#1}}\endgroup}%



% ST makrosas savo numeracijos sistemos sukurimui
\def\setvaluelist#1#2{\@tempcnta=0\relax
  \@for\@curr@val:=#2\do{%
     \advance\@tempcnta by1\relax
     \expandafter\protected@xdef\csname #1@item@\the\@tempcnta\endcsname{\@curr@val}%
     }%
     \expandafter\protected@xdef\csname #1@item@0\endcsname{\the\@tempcnta}%
}
\xdef\getitemvalue#1#2{\noexpand\csname #1@item@#2\endcsname}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \ead, \printead

\RequirePackage{keyval}

\def\email@text{}
\def\url@text{http://}
\def\ead@sep{;~}

% naudojame keyval paketa
\define@key{ead}{email}[true]{\def\ead@type{email}}
\define@key{ead}{url}[true]{\def\ead@type{url}}
\define@key{ead}{label}{\def\ead@label{#1}}

\DeclareRobustCommand\ead[2][label= ,email]{{%
  \def\ead@type{email}% default
  \setkeys{ead}{#1}%
  \def\texttildelow{\noexpand\texttildelow}%
   \protected@edef\@currentlabel{#2}%
   \safelabel{\ead@label @\ead@type}}}

\newif\ifnot@ead@star

\DeclareRobustCommand{\printead}{\@ifstar{\not@ead@starfalse\@printead}{\not@ead@startrue\@printead}}

\def\@printead{\@ifnextchar[{\@tempswatrue\@@printead}{\@tempswafalse\@@printead[]}}

\def\@@printead[#1]#2{%
   \if@tempswa% []
     {\ead@size #1}%
   \else% 
     \def\ead@type{email}%
     \def\ead@prefix{mailto:}%
     \let\ead@text\email@text%
     \let\@ead@sep\relax%
     \@for\ead@ref:=#2\do{%
         \@ead@sep\let\@ead@sep\ead@sep%
         \@ifundefined{r@\ead@ref @url}{}{\let\ead@text\url@text\def\ead@type{url}\def\ead@prefix{http://}}%
         \ifnot@ead@star\ead@text\fi{\ead@size\def\null{}\ims@href{\ead@prefix\saferef{\ead@ref @\ead@type}}{\saferef{\ead@ref @\ead@type}}}%
         \let\ead@text\relax}%                        
   \fi%
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \copyrightowner, \corref

\def\copyrightowner#1{\def\@copyrightowner{#1}}


% for corresponding author
\def\corref#1{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normaltext, nohyphen, today

% normalus tekstas (justify)
\def\normaltext{\let\\=\@normalcr%
  \leftskip\z@ \@rightskip\z@ \rightskip\@rightskip%
  \parfillskip\@flushglue}

% skiemenavimo isjungimas
\def\nohyphen{\pretolerance=\@M \tolerance=\@M \hyphenpenalty=\@M \exhyphenpenalty=\@M}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\pagestyle{ectaheadings} \thispagestyle{copyright}


\if@twoside

  \def\ps@ectaheadings{%
        \def\@oddfoot{\hfill\info@line}%
        \let\@evenfoot\@oddfoot%
        \def\@evenhead{\runninghead@size\rlap{\pagenumber@size\thepage}\hfill\leftmark\hfill}%
        \def\@oddhead{\runninghead@size\hfill\rightmark\hfill\llap{\pagenumber@size\thepage}}}
\else
   \def\ps@ectaheadings{%
         \def\@oddfoot{\hfill\info@line}%
         \let\@evenfoot\@oddfoot%
         \def\@evenhead{\runninghead@size\hfill\leftmark/\rightmark\hfill\llap{\pagenumber@size\thepage}}%
         \def\@oddhead{\runninghead@size\hfill\leftmark/\rightmark\hfill\llap{\pagenumber@size\thepage}}}
\fi

\def\ps@copyright{\let\@mkboth\@gobbletwo%
  \def\@evenhead{\parbox[t]{\textwidth}{\copyright@size\copyright@text}}%
  \let\@oddhead\@evenhead%
  \def\@oddfoot{\hfill\pagenumber@size\thepage\hfill\llap{\info@line}}%
  \let\@evenfoot\@oddfoot}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% url@fmt


\def\journalurl#1{\def\journal@url{#1}}
\def\paperurl#1{\def\paper@url{#1}}

\def\doi#1{%
  \gdef\@doi{#1}%
  \gdef\doi@text{\url@fmt{DOI: }{\ttfamily}{#1}{\doi@base\@doi}}%
}

\def\@doi{}

\def\doi@base{http://dx.doi.org/}

% {url}{text}
\def\ims@href#1#2{#2}

% {prefix}{font}{text}{url}

\def\url@fmt#1#2#3#4{%
   \edef\@tempa{#3}%
   \ifx\@tempa\@empty%
   \else%
     #1{#2\ims@href{#4}{#3}}%
   \fi}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LIST ENVIRONMENTS

\parsep\z@
\topsep\z@
\partopsep\z@
\itemsep\z@
\labelsep.5em

\def\@listI{\leftmargin\leftmargini
            \parsep\z@
            \topsep\z@
            \itemsep\z@}

\def\list@parindent{1pc}

% quotation
\let\quotation@size\footnotesize
\def\quotation@itemindent{\list@parindent}
\def\quotation@parindent{\list@parindent}
\def\quotation@leftmargin{\list@parindent}
\def\quotation@rightmargin{\list@parindent}
\let\quotation@topsep\smallskipamount

\def\quotation{%
        \list{}{\quotation@size%
        \listparindent\quotation@parindent%
        \itemindent   \quotation@itemindent%
        \rightmargin\quotation@rightmargin   \leftmargin\quotation@leftmargin%
        \partopsep\z@ \topsep\quotation@topsep \parsep\z@%
                        }%
        \item[\Q@strut]\relax}

\def\endquotation{\endlist}

\def\Q@strut{\leavevmode\hbox{\vrule height9pt depth1pt width0pt}}

% quote
\let\quote@size\footnotesize
\def\quote@indent{\z@}
\def\quote@leftmargin{2pc}
\def\quote@rightmargin{2pc}
\let\quote@topsep\smallskipamount

\def\quote{%
        \list{}{\quote@size%
        \listparindent\quote@indent%
        \itemindent \listparindent%
        \rightmargin\quote@rightmargin   \leftmargin\quote@leftmargin%
        \partopsep\z@ \topsep\quote@topsep \parsep\z@%
                       }%
        \item\relax}

\def\endquote{\endlist}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% table, figure

% table

\def\thetable{\@Roman\c@table}
\def\fnum@table{\tablename~\thetable}
\setlength\belowcaptionskip{4\p@}

\renewenvironment{table}
               {\let\@makecaption\@maketablecaption\@float{table}}
               {\end@float}
\renewenvironment{table*}
               {\let\@makecaption\@maketablecaption\@dblfloat{table}}
               {\end@dblfloat}

\long\def\@maketablecaption#1#2{%
    \parbox[t]{\hsize}{
      \centering\footnotesize\scshape
      \MakeUppercase{#1}\endgraf\vskip3\p@
      #2\par}
  \vskip\belowcaptionskip}

% figure : use \@makecaption:
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{\hskip\true@parindent\textsc{#1}.--- #2}%
  \ifdim \wd\@tempboxa >\hsize
    \hskip\true@parindent\textsc{#1}.--- #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\hskip-\true@parindent\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}

%
\def\@floatboxreset{%
        \reset@font
        \footnotesize
        \@setminipage
        \centering
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOOTNOTE

\renewcommand\@makefntext[1]{%
    \parindent\true@parindent\@makefnmark #1}

\def\@makefnmark{\@textsuperscript{\normalfont\@thefnmark}}%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FRONTMATTER 
% FRONT MATTER FORMATTING PARAMETERS

\setattribute{frontmatter} {style} {\centering}
\setattribute{title}       {style} {\centering}
\setattribute{author}      {style} {\centering}
\setattribute{address}     {style} {}
\setattribute{affiliation} {style} {\normaltext}
\setattribute{abstract}    {style} {\normaltext}
\setattribute{keyword}     {style} {\normaltext}
\setattribute{history}     {style} {\centering}

% FRONT MATTER SKIPS
\setattribute{title}       {skip} {-6\p@}
\setattribute{authors}     {skip} {\medskipamount}
\setattribute{dedicated}   {skip} {\medskipamount}
\setattribute{address}     {skip} {\smallskipamount}
\setattribute{affiliation} {skip} {\medskipamount}
\setattribute{abstract}    {skip} {\true@parindent}
\setattribute{keyword}     {skip} {\true@parindent}
\setattribute{history}     {skip} {\true@parindent}
\setattribute{frontmatter} {cmd}  {\vskip\bigskipamount
                                   \@afterindentfalse\@afterheading}
% FRONT MATTER DIMENSIONS
\if@xiisize
  \setattribute{abstract}   {width} {29pc}
  \setattribute{keyword}    {width} {29pc}
\else
  \setattribute{abstract}   {width} {26pc}
  \setattribute{keyword}    {width} {26pc}
\fi
\setattribute{abstract}  {indent} {\true@parindent} 
\setattribute{keyword}   {indent} {\true@parindent} 
\setattribute{history}    {width} {\textwidth}

% FRONT MATTER FONTS 
\setattribute{dochead}    {size} {\normalsize}
\setattribute{title}      {size} {\fontsize{11}{13}\selectfont}
\setattribute{author}     {size} {\normalsize\scshape}
\setattribute{address}    {size} {\normalsize\itshape}
\setattribute{affiliation}{size} {\normalsize\itshape}
\setattribute{dedicated}  {size} {\normalsize\itshape}
\setattribute{ead}        {size} {\upshape\ttfamily}
\setattribute{abstract}   {size} {\footnotesize\upshape\mdseries}
\setattribute{abstractname}{size} {}
\setattribute{keyword}    {size} {\footnotesize\upshape\mdseries}
\setattribute{keywordname}{size} {\scshape}
\setattribute{history}    {size} {\footnotesize\itshape}
\setattribute{copyright}  {size} {\footnotesize\centering}
\setattribute{runninghead}{size} {\footnotesize}
\setattribute{pagenumber} {size} {\normalsize\upshape}
\setattribute{appendix}   {size} {\footnotesize}
\setattribute{thebibliography}{size}{\footnotesize}

% FRONT MATTER CASE
\setattribute{dochead}    {case} {\MakeUppercase}
\setattribute{title}      {case} {\MakeUppercase}
\setattribute{runninghead}{case} {\MakeUppercase}

% TEXT, etc.
\setattribute{volume}     {text} {Vol.}
\setattribute{doi}        {text} {\url@fmt{url: }{\ttfamily}{\paper@url}{\paper@url}}
\setattribute{copyright}  {text} {\url@fmt{}{\itshape}{\journal@name}{\journal@url}}
\setattribute{infoline}   {text} {\lower12pt \hbox{\scriptsize\ttfamily\ectafmt@name\ ver. \ectafmt@version\ file: \jobname.tex\ date: \today}}
\setattribute{copyright} {owner} {$\copyright$~\@copyrightyear \@copyrightowner}
\setattribute{author}   {prefix} {}
\setattribute{keyword} {postfix} {\unskip.}

\def\abstractname{}

% HISTORY
\setattribute{history}  {prefix}  {}
\setattribute{history}  {postfix} {.}
\setattribute{received} {prefix}  {Manuscript received~}
\setattribute{received} {postfix} {}
\setattribute{revised}  {prefix}  {; final revision received~}
\setattribute{revised}  {postfix} {}
\setattribute{accepted} {prefix}  {; accepted~}
\setattribute{accepted} {postfix} {}


% COUNTERS, ETC
\newcounter{author}
\newcounter{address}
\newdimen\sv@mathsurround
\def\author@num{0}

% RUNNING HEAD
\def\runtitle#1{\gdef\@runtitle{\runninghead@case{#1}}}                      
\def\runauthor#1{{\def\etal{et al.}\gdef\@runauthor{\runninghead@case{#1}}}} 

\newdimen\sv@parindent
\sv@parindent\parindent

\newbox\fm@box
\newdimen\fm@size

\let\hy@frontmatter\relax
\let\hy@endfrontmatter\relax

\newif\if@frontmatter  \@frontmatterfalse

\def\frontmatter{%
  \global\c@author\z@
  \global\c@address\z@
  \@frontmattertrue
  \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
  \sv@mathsurround\mathsurround \m@th
  \parindent\z@
  \hy@frontmatter
  \global\let\maketitle\relax
  \open@fm \ignorespaces}

% ENDFRONTMATTER

\def\endfrontmatter{%
  \global\@topnum\z@
  \set@pagerange
  \@ifundefined{@runtitle}{\def\@runtitle{\runninghead@case{\safe@title}}}{}
  \@ifundefined{@runauthor}{\def\@runauthor{\runninghead@case{\the\authors@list}}}{}
  \if@ectadraft
    \markboth{}{}%
  \else
    \markboth{\@runauthor}{\@runtitle}%
  \fi
  \thispagestyle{copyright}%
%
  \put@fmt@data%
%
  \close@fm
  \write@pdfinfo{\hy@author}{\the\authors@list}
  \write@pdfinfo{\hy@keywords}{\the\keywords@list}
%
  \write\@mainaux{\string\gdef\string\author@num{\the\c@author}}
%
  \hy@endfrontmatter
  \global\mathsurround\sv@mathsurround
  \global\c@footnote\z@
  \global\let\@thanks\@empty  
  \global\let\title\relax       
  \global\let\author\relax
  \aftergroup\frontmatter@cmd
  }

\def\put@fmt@data{%
  \copyright@fmt%
  \@thanks%
  \abstract@fmt%
  \keyword@fmt}

\newdimen\t@xtheight
\def\init@settings{
\splittopskip=\topskip \splitmaxdepth=\maxdepth
\t@xtheight\textheight \advance\t@xtheight-\splittopskip}

\def\no@harm{\let\thanks=\@gobble\let\thanksref=\@gobble\let~\space\def\ead[##1]##2{}\let\\=\@empty \def\protect{\noexpand\protect\noexpand}}

\def\open@fm{%
  \global\setbox\fm@box=\vbox\bgroup
  \hsize=\textwidth
  \frontmatter@style}
%
\def\close@fm{%
  \par \egroup
  \print@titlepage
  \fm@size=\dp\fm@box \advance\fm@size by \ht\fm@box
  \@whiledim\fm@size>\t@xtheight \do{%
    \global\setbox\@tempboxa=\vsplit\fm@box to \t@xtheight
    \unvbox\@tempboxa %\newpage
    \fm@size=\dp\fm@box \advance\fm@size by \ht\fm@box}
  \if@twocolumn
    \emergencystretch=1pc \twocolumn[\unvbox\fm@box]
  \else
    \unvbox\fm@box
  \fi}

% DOCHEAD
\def\dochead#1{%
  \bgroup
    \dochead@size
    \leavevmode\vphantom{\strut}\dochead@case{#1}\par
  \egroup
  \setattribute{title}{skip}{8\p@}}

% TITLE

\def\title#1{%
  \vglue\title@skip%
% check if we are in {frontmatter}
  \if@frontmatter \else\hy@frontmatter\fi
  \bgroup%
    \no@harm%
    \let\protect\relax%
    \xdef\@argi{#1}%
    \xdef\safe@title{#1}%
  \egroup%
  \write@pdfinfo{\hy@title}{\@argi}
  \if@frontmatter%
    \bgroup%
      \title@style\title@size\title@case{#1}\par%
    \egroup%
  \else%
    \gdef\@title{#1}%
  \fi}

% AUG - author block
\def\smart@par{\ifhmode\par\fi}
\newenvironment{aug}{}{\smart@par}

% AUTHOR

\def\and{\unskip~and~}

\def\author{\@ifnextchar[{\author@fmt}{\author@fmt[]}}

\def\author@fmt[#1]#2{%
  \stepcounter{author}%
  \author@fmt@init%
  \let\author@fmt@init\relax%
  \bgroup% 
     \def\degs##1{##1}\def\fnms##1{##1}\def\inits##1{##1}\def\snm##1{##1}\def\roles##1{##1}%
     \@tempcnta=\author@num\relax%
     \ifnum\c@author=\@tempcnta \def\author@sep{ and }\else \def\author@sep{, }\fi%
     \ifnum\c@author=1\def\author@sep{}\fi%
     \addto@authors@list{\author@sep #2}%
     \if@frontmatter%
       \author@sep #2\thanksref{#1}%
     \else \gdef\@author{#2}\fi%
  \egroup%
  \ignorespaces}

\def\author@fmt@init{%
   \vskip\authors@skip%
   \leavevmode\author@style\author@size\author@prefix}

\let\author@fmt@init@def\author@fmt@init

% DEDICATED
\def\dedicated#1{%
  \vskip\dedicated@skip
  \bgroup
    \dedicated@size #1\par
  \egroup}

% ADDRESS

\def\address{\@ifnextchar[{\address@fmt}{\address@fmt[default]}}
%
% \def\address@ref@A=1
% \newtoks\@address@1
% \@address@1={#2}

\def\address@fmt[#1]#2{%
  \stepcounter{address}%
  \expandafter\protected@xdef\csname address@ref@#1\endcsname{\the\c@address}%
  \expandafter\newtoks\csname @address@\the\c@address\endcsname
  \expandafter\expandafter\expandafter\global\csname @address@\the\c@address\endcsname={#2}}


\let\safe@phantomsection\relax
% print address by number: \printaddressnum{1}
\def\printaddressnum#1{%
\xdef\@tmp{#1}%
\bgroup
\@ifundefined{@address@#1}{\@latex@error{Klaida: nera adreso su numeriu '#1'!}{??}}{
  \address@size
  \ifnum#1=1%
    \leavevmode\safe@phantomsection\addcontentsline{toc}{section}{Author's addresses}%
  \fi%
  \expandafter\expandafter\expandafter\the\csname @address@\@tmp\endcsname
  \par}
\egroup
}

%
\def\print@address@and{{\centering\itshape and\par}}

% print all addresses:
\def\address@par{\par\vskip3pt}

\def\printaddresses{%
\vskip\address@skip%
\ifnum\c@address=1\let\print@address@and\relax\fi%
\@tempcnta=0%
\bgroup%
\@whilenum{\@tempcnta<\c@address}%
  \do{%
    \advance\@tempcnta\@ne
    \ifnum\@tempcnta=\c@address \print@address@and\fi %paskutinis narys
    \expandafter\printaddressnum{\the\@tempcnta}%
     }
\egroup
\global\let\printaddresses\relax
}


%%%%%%%%%%% new layout adrress as footnote:

\if@addressasfootnote

  \newcounter{addressref}
  \def\theaddressref{{\upshape\alph{addressref}}}

  \long\def\safe@footnotetext#1{\insert\footins{%
      \reset@font\footnotesize
      \interlinepenalty\interfootnotelinepenalty
      \splittopskip\footnotesep
      \splitmaxdepth \dp\strutbox \floatingpenalty \@MM
      \hsize\columnwidth \@parboxrestore
      \color@begingroup
        \def\@thefnmark{}%
        \@makefntext{%
          \rule\z@\footnotesep\ignorespaces #1\@finalstrut\strutbox}%
      \color@endgroup}}%

  \def\address{\@ifnextchar[{\address@fmt}{\address@fmt[]}}

  \def\address@fmt[#1]#2{%
     \ifx.#1.%
       \safe@footnotetext{#2}
     \else%
       \usethankscounter{addressref}%
       \bgroup
         \def\\{\hfill\break}
         \thankstext{#1}{#2}%
       \egroup
       \restorethankscounter%
     \fi%
                        }

  \let\printaddresses\relax

\fi




% AFFILIATION
\def\affiliation#1{%
  \smart@par%
  \let\author@fmt@init\author@fmt@init@def
  \vskip\affiliation@skip
  \bgroup
    \affiliation@size #1\par
  \egroup}


% COPYRIGHTOWNER

\def\copyright@fmt{%
  \@ifundefined{@copyrightowner}{}{%
  \bgroup%
    \def\thempfn{}%
    \footnotetext{\copyright@owner}%
  \egroup}}

% ABSTRACT

\newbox\abstract@box

\def\abstract{%
     \global\setbox\abstract@box=\vbox\bgroup%
     \@tempdima\textwidth%
     \advance\@tempdima by-\abstract@width%
     \divide\@tempdima by2%
     \abstract@style%
     \leftskip\@tempdima\rightskip\@tempdima%
     \abstract@size%
     \parindent\sv@parindent%
     \noindent\hskip\abstract@indent{\abstractname@size\abstractname}\ignorespaces}
\def\endabstract{\par\egroup}

\def\abstract@fmt{%
  \ifvoid\abstract@box\else
    \vskip\abstract@skip%
    \unvbox\abstract@box
  \fi}


% HISTORY: received, revised, accepted

\def\history@exist{0}

\def\received#1{\def\@tempa{#1}\ifx\@tempa\@empty\else\gdef\@received{#1}\gdef\history@exist{1}\fi}
  \def\@received{\@nil}
\def\revised#1{\def\@tempa{#1}\ifx\@tempa\@empty\else\gdef\@revised{#1}\gdef\history@exist{1}\fi}
  \def\@revised{\@nil}
\def\accepted#1{\def\@tempa{#1}\ifx\@tempa\@empty\else\gdef\@accepted{#1}\gdef\history@exist{1}\fi}
  \def\@accepted{\@nil}

\def\empty@data{\@nil}

\def\history@fmt{%
  \ifcase\history@exist\else%
  \bgroup
    \nobreak%
    \vskip\history@skip%
    \nobreak%
    \history@style%
    \history@size%
    \leavevmode
    \history@prefix
    \ifx\@received\empty@data \else
      \received@prefix\@received \received@postfix%
    \fi
    \ifx\@revised\empty@data \else
      \revised@prefix\@revised \revised@postfix%
    \fi
    \ifx\@accepted\empty@data \else
      \accepted@prefix\@accepted \accepted@postfix%
    \fi
  \history@postfix\par%
  \egroup%
  \gdef\history@exist{0}
\fi
}

\let\printhistory\history@fmt

\def\sday#1{#1}
\def\smonth#1{\@ifundefined{month@item@#1}%
     {\@latex@error{Nera tokio menesio, kurio numeris #1!}{??}}%
     {\getitemvalue{month}{#1}}%
  }%
\def\syear#1{#1}
\setvaluelist{month}{January,February,March,April,May,June,July,August,September,October,November,December}

% KEYWORDS

\newbox\keyword@box
\newdimen\pre@kwd@depth

\setattribute{keyword}{AMS}{AMS 2000 subject classifications:}
\setattribute{keyword}{KWD}{Keywords:}
% raktas=class
\setattribute[default]{keyword}{class}{KWD}

%  \keyword@class-> KWD
%  \keyword@KWD  -> AMS 2000... 

\gdef\keyword{\@ifnextchar[{\@keyword}{\@keyword[class=KWD]}}

\gdef\@keyword[#1]{%
  \do@option@list{keyword}{#1}%
  \def\keyword@name{\csname keyword@\keyword@class\endcsname}%
  \let\kwd@sep\relax
%
  \global\setbox\keyword@box=\vbox\bgroup%
     \ifvoid\keyword@box\else%
        \unvbox\keyword@box
        \vskip-\pre@kwd@depth\vtop to\pre@kwd@depth{}%
     \fi
     \@tempdima\textwidth%
     \advance\@tempdima by-\keyword@width%
     \divide\@tempdima by2%
     \keyword@style%
     \leftskip\@tempdima\rightskip\@tempdima%
     \keyword@size%
     \parindent\sv@parindent%
     \noindent\hskip\keyword@indent{\keywordname@size\keyword@name}\space\hskip.1pt}

\def\endkeyword{\keyword@postfix\par\global\pre@kwd@depth\prevdepth\egroup}

\def\keyword@fmt{%
  \ifvoid\keyword@box\else
    \vskip\keyword@skip%
    \unvbox\keyword@box
  \fi}


% \kwd[; ]{foo}
  \def\sep{\unskip\string, }%

  \DeclareRobustCommand*\kwd{\@ifnextchar[\@kwd{\@kwd[\kwd@sep]}}%
  \def\@kwd[#1]#2{\unskip#1{#2}\addto@keywords@list{#2, }\let\kwd@sep\sep}%  


%%%%% Tittle (Cover) page

\def\print@titlepage{
  \if@titlepage
    \thispagestyle{empty}
    \bgroup
      \centering
      \renewcommand{\baselinestretch}{1.0}\large\normalsize
      {\large\journal@name\par}
      \vspace*{1pt}
      \vskip48\p@
      {\Huge\safe@title\par}
      \bigskip
      {\Large\the\authors@list\par}
      \bigskip
      \let\address@size\large
      \printcoveraddresses
      \vfill
      {\large\@date\par}
    \egroup
    \clearpage
    \thispagestyle{copyright}
    \addtocounter{page}{-1}
  \fi}

\def\printcoveraddresses{%
\ifnum\c@address=1 \let\print@address@and\relax\fi%
\@tempcnta=0%
\bgroup\parindent\z@
\@whilenum{\@tempcnta<\c@address}%
  \do{%
    \advance\@tempcnta\@ne
    \ifnum\@tempcnta=\c@address \print@address@and\fi %paskutinis narys
    \expandafter\printaddressnum{\the\@tempcnta}%
     }
\egroup
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \maketitle
% if \frontmatter is not used, we will redefine \maketitle

\def\local@maketitle{%
  \global\@topnum\z@
  \set@pagerange
%
  \@ifundefined{@runtitle}{\def\@runtitle{\runninghead@case{\safe@title}}}{}
  \@ifundefined{@runauthor}{\def\@runauthor{\runninghead@case{\the\authors@list}}}{}
%
  \markboth{\@runauthor}{\@runtitle}%
  \thispagestyle{copyright}%
%
  \print@titlepage
% 
  \bgroup%
     \title@style\title@size\title@case{\@title}\par%
  \egroup%
%
  \vskip\authors@skip%
  \bgroup    
   \leavevmode\author@style\author@size\author@prefix
   \@author\par
  \egroup
%
  \put@fmt@data%  
%
  \write@pdfinfo{\hy@author}{\the\authors@list}
  \write@pdfinfo{\hy@keywords}{\the\keywords@list}
  \hy@endfrontmatter
  \global\mathsurround\sv@mathsurround
  \global\c@footnote\z@
  \global\let\@thanks\@empty  
  \let\title\@gobble
  \let\author\@gobble
  \let\frontmatter\relax \let\endfrontmatter\relax
  \let\@maketitle\relax  \let\@@maketitle\relax
  \normalfont\normaltext
  \frontmatter@cmd
%
  \ifvoid\abstract@box\g@addto@macro\endabstract\abstract@fmt\fi
  \ifvoid\keyword@box\g@addto@macro\endkeyword\keyword@fmt\fi
  }

\AtBeginDocument{\let\maketitle\local@maketitle}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% interaction with hyperref

\def\test@hyperref{\@ifundefined{Hy@SetCatcodes}{}{\ectaart@hyperref@settings}}

\AtBeginDocument{\test@hyperref}

\def\ectaart@hyperref@settings{%
% from hyperref maketitle
  \def\hy@frontmatter{%
    \let\H@@origfootnotemark\@footnotemark
    \let\H@@origfootnotetext\@footnotetext
    \let\@footnotemark\H@@footnotemark
    \let\@footnotetext\H@@footnotetext}
  \def\hy@endfrontmatter{
    \ifx\@footnotemark\H@@footnotemark
      \let\@footnotemark\H@@origfootnotemark
    \fi
    \ifx\@footnotetext\H@@footnotetext
      \let\@footnotetext\H@@origfootnotetext
    \fi}
% hooks for the \thanksref, \thankstext:
  \def\thankref@hyperlink##1{%
     \edef\@tempx{##1thanks}%
     \hbox{\hyperlink{##1}{\saferef{\@tempx}}}}
  \def\thanks@hypertarget##1{\smash{\raise\baselineskip\hbox{\protect\hypertarget{##1}{}}}}
% redefine pagenumbering
  \let\pagenumbering\ims@pagenumbering
% activate href
   \let\ims@href\href%
   \let\safe@phantomsection\phantomsection
% put document info
   \def\write@pdfinfo##1##2{\protected@write\@auxout{\no@harm}{\string\gdef\string##1{##2}}}
   \@ifundefined{hy@title}{}{\pdfstringdef\@pdftitle{\hy@title}}
   \@ifundefined{hy@author}{}{\pdfstringdef\@pdfauthor{\hy@author}}
   \@ifundefined{hy@keywords}{}{\pdfstringdef\@pdfkeywords{\hy@keywords}}
% MathSciNet:
   \def\MR##1{\href{http://www.ams.org/mathscinet-getitem?mr=##1}{MR##1}}
}

\def\write@pdfinfo#1#2{}

\newtoks\authors@list
\def\addto@authors@list#1{%
  \begingroup%
    \no@harm%
    \xdef\@act{\global\noexpand\authors@list{\the\authors@list#1}}\@act%
  \endgroup}

\newtoks\keywords@list
\def\addto@keywords@list#1{%
  \begingroup%
    \no@harm%
    \xdef\@act{\global\noexpand\keywords@list{\the\keywords@list#1}}\@act%
  \endgroup}


% MathSciNet:
\def\MR#1{MR#1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% section, subsection etc.
% we do not like article appearance:

\let\section@cnt@size\normalsize

% paskutine eilute centruojame *lapo* atzvilgiu, o ne pirmos
% eilutes teksto atzvilgiu

\def\fl@hangfrom#1{\noindent {#1}}
\def\TO@fl{\let\@hangfrom\fl@hangfrom}


\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {\bigskipamount}%
                                   {\smallskipamount}%
                                   {\footnotesize\baselineskip12\p@\centering\TO@fl\MakeUppercase}}

\renewcommand\subsection{\@startsection {subsection}{2}{\z@}%
                                   {\bigskipamount}%
                                   {\smallskipamount}%
                                   {\itshape\centering\TO@fl}}

\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                   {\medskipamount}%
                                   {\smallskipamount}%
                                   {\itshape\raggedright}}

\renewcommand\paragraph{\@startsection{paragraph}{4}{\parindent}%
                                   {\medskipamount}%
                                   {\true@parindent}%
                                   {\itshape}}

\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                   {0.1pt}%
                                   {-1em}%
                                   {\itshape}}



% Format for the counter:
  \def\section@numbersep{.}
  \def\subsection@numbersep{.}
  \def\subsubsection@numbersep{.}
  \def\paragraph@numbersep{.}
  \def\subparagraph@numbersep{.}

  \def\@seccntformat#1{{\section@cnt@size\upshape\csname #1@prefix\endcsname\csname the#1\endcsname\csname#1@numbersep\endcsname\enspace}}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% appendix


  \renewcommand\appendix{\par
    \printaddresses%
    \history@fmt%
    \let\section@cnt@size\appendix@size%
    \let\old@section\section%
    \def\section{\@ifnextchar*{\@appsectionstar}{\@appsectionnostar}}%
    \def\section@prefix{\appendixname\ }%
    \def\section@numbersep{:}%
    \setcounter{section}{0}%
    \setcounter{subsection}{0}%
    \gdef\thesection{\@Alph\c@section}
    \appendix@size
}

  \def\@appsectionstar*#1{%
     \old@section*{#1}%
     \setcounter{section}{1}%
     \addcontentsline{toc}{section}{#1}
}

  \def\@appsectionnostar#1{%
     \ifx.#1.% 
       \def\section@numbersep{}\old@section[\appendixname\ \thesection]{}%
     \else%
       \def\section@numbersep{:}\old@section{#1}%
     \fi}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% thebibliography = article

\def\ectafix@thebibliography{%
  \let\old@thebibliography\thebibliography
  \def\thebibliography{%
  \printaddresses%
  \history@fmt%
  \old@thebibliography}}

\AtBeginDocument{\ectafix@thebibliography}

\def\bysame{\@ifnextchar\bgroup{\@bysame}{\@bysame{}}}
\def\@bysame#1{\vrule height 1.5pt depth -1pt width 3em \hskip 0.5em\relax}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fix if natbib is loaded:

\def\ectafix@natbib{%
   \let\bibfont\thebibliography@size
   \setlength\bibsep{0pt}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fix if amsthm isloaded:

\def\ectafix@amsthm{%
     \newtheoremstyle{plain}     {\medskipamount}{\medskipamount}{\itshape}{}{\scshape}{}{1em}{}
     \newtheoremstyle{definition}{\medskipamount}{\medskipamount}{\normalfont}{}{\scshape}{}{1em}{}
     \newtheoremstyle{remark}    {\medskipamount}{\medskipamount}{\normalfont}{}{\scshape}{}{1em}{}
%
   \def\qed{\hfill\hbox{\textit{Q.E.D.}}}
%
\renewenvironment{proof}[1][\proofname]{\par
  \pushQED{\qed}%
   \setlength\labelsep{1em}
  \normalfont \topsep6\p@\@plus6\p@\relax
  \trivlist
  \item[\hskip\labelsep
        \scshape
    ##1\@addpunct{:}]\ignorespaces
}{\popQED\endtrivlist\@endpefalse}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check for the loaded classes:

\AtBeginDocument{%
   \@ifpackageloaded{natbib}{\ectafix@natbib}{}%
   \@ifpackageloaded{amsthm}{\ectafix@amsthm}{}%
   \@ifpackageloaded{amsmath}{\@mathmargin=42\p@ minus42\p@}{}%
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% number lines

\ifnumberlines@

   \long\def\l@addto@macro#1#2{%
       \toks@\expandafter{#1#2}%
       \edef#1{\the\toks@}}

   \setattribute{numberlines}{size}{\scriptsize\ttfamily}
   \setattribute{numberlines}{skip}{64\p@}

   \def\numberlines@hook{%
       \l@addto@macro\@evenhead\put@numberlines@box%
       \l@addto@macro\@oddhead\put@numberlines@box}

   \g@addto@macro\ps@ectaheadings\numberlines@hook
   \g@addto@macro\ps@copyright\numberlines@hook

   \newbox\numberlines@box
   \newskip\numberlines@box@skip

   \def\set@numberlines@box{%
     \setlength\numberlines@box@skip\headsep
     \addtolength\numberlines@box@skip{5\p@}
   %
     \setbox\numberlines@box\vtop to\textheight{%
       \parindent\z@    
       \vskip\z@   
       \@tempcnta=0
       \@tempdima=\z@
       \loop
         \advance\@tempcnta by1
         \advance\@tempdima by\baselineskip
         \hbox to\textwidth{%
            \llap{\numberlines@size\the\@tempcnta\kern\numberlines@skip}
            \hfill
            \rlap{\numberlines@size\kern\numberlines@skip\the\@tempcnta}}
       \ifdim\@tempdima<\textheight\repeat
       \vss
     }%
   %
       \ht\numberlines@box\z@
       \dp\numberlines@box\z@
   }

   \def\put@numberlines@box{\lower\numberlines@box@skip\hbox to\z@{\hss\copy\numberlines@box}}

   \AtBeginDocument{\set@numberlines@box}

\fi

%**************** INICIALIZATION

% print addresses if they were not printed:
\AtEndDocument{
  \printaddresses%
  \history@fmt}

\@twosidetrue
\pagenumbering{arabic}
\frenchspacing
\init@settings
\pagestyle{ectaheadings}

\endinput
%%
%% End of file `ectaart.cls'.
