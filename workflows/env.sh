#!/bin/bash

export NMR_IO_TIMEOUT=0
export NMR_IO_SELECT=0
export NMR_AUTOSWAP=1

alias u='chmod a+rx *.com *.tcl'

#set   history = 256
#alias h      'history'
#alias rm     'rm -i'
#alias cd     'cd \!*; set prompt = "$cwd% "'
#alias dirs   'ls -l | egrep -e ^d'
#alias links  'ls -l | egrep -e ^l'
prefix="/usr/local/nmrpipe"
export MANPATH="$prefix/man:/usr/share/man:${MANPATH}"
export PATH="$prefix/nmrbin.linux9:$prefix/com:$PATH"

export LD_LIBRARY_PATH="$prefix/nmrbin.linux9/lib:${LD_LIBRARY_PATH}"

export OPENWINHOME="$prefix/nmrbin.linux9/openwin"

export NMRCHECK="ALL"
export NMRBASE="$prefix"
export NMRBINTYPE="linux9"
export NMRTXT="$prefix/nmrtxt"
export NMRBIN="$prefix/nmrbin.linux9"
export TCLPATH="$prefix/com"
export TALOS_DIR="$prefix/talos"
export TALOSP_DIR="$prefix/talosplus"
export SPARTAP_DIR="$prefix/spartaplus"
export PROMEGA_DIR="$prefix/promega"


export NMR_TCLTK8="TRUE"

if [[ $NMR_TCLTK8 ]] ; then
   export TCL_LIBRARY="$prefix/nmrtcl/tcl8.4"
   export TK_LIBRARY="$prefix/nmrtcl/tk8.4"
   export BLT_LIBRARY="$prefix/nmrtcl/blt2.4"
   export NMRPIPE_TCL_LIB="$prefix/nmrtcl/tcl8.4"
   export NMRPIPE_TK_LIB="$prefix/nmrtcl/tk8.4"
   export NMRPIPE_BLT_LIB="$prefix/nmrtcl/blt2.4"
else
   export TCL_LIBRARY="$prefix/nmrtcl/tcl7.6"
   export TK_LIBRARY="$prefix/nmrtcl/tk4.2"
   export BLT_LIBRARY="$prefix/nmrtcl/blt2.4"
   export NMRPIPE_TCL_LIB="$prefix/nmrtcl/tcl7.6"
   export NMRPIPE_TK_LIB="$prefix/nmrtcl/tk4.2"
   export NMRPIPE_BLT_LIB="$prefix/nmrtcl/blt2.4"
fi

# setenv NMRPIPE_SMALLFONT  "-adobe-helvetica-medium-r-*-*-*-100-*-*-*-*-*-*"
# setenv NMRPIPE_BIGFONT    "-adobe-helvetica-medium-r-*-*-*-180-*-*-*-*-*-*"
# setenv NMRPIPE_STDFONT    "-adobe-helvetica-medium-r-*-*-*-120-*-*-*-*-*-*"
# setenv NMRPIPE_BOLDFONT   "-adobe-helvetica-bold-r-*-*-*-120-*-*-*-*-*-*"
# setenv NMRPIPE_FIXEDFONT  "-*-courier-medium-r-*-*-*-120-*-*-*-*-*-*"
# setenv NMRPIPE_TYPINGFONT "-*-courier-medium-r-*-*-*-100-*-*-*-*-*-*"
# setenv NMRPIPE_AXISFONT   "-*-lucidatypewriter-bold-r-*-*-*-120-*-*-*-*-*-*"

if [[ $NMRFONTSET ]] ; then
   exit 0
fi

if [[ ! -d $prefix/XView/fonts/gz ]] ; then
   exit 0
fi

if [[ $DISPLAY ]] ; then
   fontInfo=( `xset -q | fgrep XView/fonts | wc` )

   if [[ $fontInfo[1] == 0 ]] ; then
      xset fp+ $prefix/XView/fonts/gz/75dpi/
      xset fp+ $prefix/XView/fonts/gz/misc/
      xset fp  rehash
     echo NMR Fonts: $prefix/XView/fonts/gz
  else
     echo NMR font path is already set.
  fi
fi

export NMRFONTSET=TRUE

export  PROMPT_COMMAND='printf "\033]0;%s@%s:%s NMRPIPE\007" "${USER}" "${HOSTNAME%%.*}" "${PWD/#$HOME/~}"'
