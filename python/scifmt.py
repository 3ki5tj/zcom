#!/usr/bin/env python



''' write a number in scientific format '''



import math



class scifmt:
  def __init__(me, x, err = 0):
    me.x = x # value

    # the sign
    if x > 0: me.sgn = 1
    elif x == 0: me.sgn = 0
    else: me.sgn = -1

    # absolute value
    me.mag = me.x * me.sgn

    # exponent
    me.exp = 0
    if me.sgn != 0:
      me.exp = int(math.log10(me.mag) + 10000) - 10000

    # x = me.mag x 10^{me.exp}
    me.mag /= math.pow(10, me.exp) # significand

    # error
    if err >= 0: me.err = err
    else: me.err = -err



  def geterr(me, err, errmax):
    ''' write error as c x 10^{xp} with c <= errmax '''
    xp = int(math.log10(err) + 1000) - 1000 + 1
    while int(err/10**xp + .5) < errmax:
      xp -= 1
    xp += 1
    c = int(err/10**xp + .5)
    return c, xp



  def getxerr(me, xmag, xexp, err, errmax, expmin = -2, expmax = 2):
    ''' write x(err) x 10^{xp} '''

    # write me.err = interr * 10^xperr
    interr, xperr = me.geterr(err, errmax)

    diffxp = xexp - xperr
    mul = 10**diffxp
    intx = int(xmag * mul + .5)
    if diffxp > 0:
      nx = "%.*f" % (diffxp, 1.*intx/mul)
      nxp = xexp
    elif diffxp <= 0:
      nx = intx
      nxp = xperr
    return nx, interr, nxp


  def shiftdigits(me, xpmin, xpmax):
    # shift the decimal mark
    mag, xp = me.mag, me.exp
    if xp <= xpmax and xp >= xpmin:
      mag *= 10**xp
      xp = 0
    return mag, xp


  def text(me, fmtmag = "%s", errmax = 50, xpmin = -2, xpmax = 2):
    ''' output to text '''
    strsgn = "- +"[me.sgn + 1]

    mag, nxp = me.shiftdigits(xpmin, xpmax)

    if me.err == 0: # no error
      s = strsgn + fmtmag % mag
    else: # has error
      nx, interr, nxp = me.getxerr(mag, nxp, me.err, errmax)
      s = strsgn + "%s(%s)" % (nx, interr)
    if nxp != 0: s += "x10^%s" % nxp
    return s



  def __str__(me):
    return me.text()



  def html(me, fmtmag = "%s", errmax = 50):
    ''' output to HTML format '''
    if me.exp >= 0: myexp = str(me.exp)
    else: myexp = "&minus;" + str(me.exp)
    strsgn = ["&minus;", "&nbsp;", "+"][me.sgn + 1]
    if me.err == 0:
      nxp = me.exp
      s = strsgn + fmtmag % me.mag
    else:
      nx, interr, nxp = me.getxerr(me.mag, me.exp, me.err, errmax)
      s = strsgn + "%s(%s)" % (nx, interr)
    if nxp != 0: s += "&times;10<sup>%s</sup>" % nxp
    return s



  def latex(me, fmtmag = "%s", errmax = 50):
    ''' output to Latex '''
    strsgn = ["-", "", "+"][me.sgn + 1]
    if me.err == 0:
      nxp = me.exp
      s = "$%s%s" % (strsgn, me.mag)
    else:
      nx, interr, nxp = me.getxerr(me.mag, me.exp, me.err, errmax)
      s = "$%s%s(%s)" % (strsgn, nx, interr)
    if nxp != 0: s += "\\times 10^{%s}$" % nxp
    else: s += "$"
    return s



if __name__ == "__main__":
  arr = [scifmt(132.2312), scifmt(0), scifmt(8.9), scifmt(-0.012),
         scifmt(12.345, 0.6), scifmt(0.6, 1.3), scifmt(3.2, 0.005),
         scifmt(9.832, 0.012), scifmt(3, 2508)]
  for x in arr:
    print x, "\t", x.html(), "\t", x.latex()
  print scifmt(0.3, 250.8).text(errmax = 3999)
