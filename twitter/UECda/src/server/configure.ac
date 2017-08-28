#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT(tndhms, 0.09, daihinmin@ice.uec.ac.jp)
AM_INIT_AUTOMAKE(tndhms, 0.09, daihinmin@ice.uec.ac.jp)
AC_CONFIG_SRCDIR([src/define.h])
AC_CONFIG_HEADER([src/config.h])

# Checks for programs.
AC_PROG_CC
# AC_PATH_X
AC_PATH_XTRA

# Checks for libraries.
AC_DEFUN([LIB_SOCKET_NSL],
[
        AC_SEARCH_LIBS([gethostbyname], [nsl])
        AC_SEARCH_LIBS([socket], [socket], [], [
                AC_CHECK_LIB([socket], [socket], [LIBS="-lsocket -lnsl $LIBS"],
                [], [-lnsl])])
])
LIB_SOCKET_NSL
LDFLAGS="$LDFLAGS $X_LIBS"
AC_CHECK_LIB([X11], [XmbDrawString])
AC_CHECK_LIB([Xpm], [XpmCreatePixmapFromData])
#AC_CHECK_LIB([Xt], [main],,,$X_LIBS $LIBS)
#AC_CHECK_LIB([Xaw], [main],,,$X_LIBS $LIBS)
#AC_CHECK_LIB([], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([arpa/inet.h inttypes.h netdb.h netinet/in.h stdlib.h string.h strings.h sys/param.h sys/socket.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([bzero inet_ntoa memset socket XCopyArea XpmReadFileToPixmap],,AC_MSG_ERROR(not found.))

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
