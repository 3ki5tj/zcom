#!/usr/bin/env python
import os, sys, re

''' script to convert .def to input '''

file = "at_t.def"
var_prefix = ""
key_prefix = ""
flag_prefix = ""
xm_prefix = "AT_"
cond0_common = ""
def sub(line):
  line = re.sub(r"AT_SET_TRACE", r'wtrace_buf("%@filename=%s", AT_(trace_file)) == 0', line)
  line = re.sub(r"AT_BOOST_FORMULA",  r"H = H0 + [shift0 + alpha0 * (T/Tref -1)] * H1", line)
  line = re.sub(r"AT_BOOST_MODE", r"@boost_mode", line)
  return line

r'''
file = "mb_t.def"
var_prefix = ""
key_prefix = ""
flag_prefix = "MBEST_"
xm_prefix = "MB_"
cond0_common = ""
def sub(line):
  line = re.sub(r"MB_MBDEL_VALID\((.*?)\)", r"@mbin_delta > @beta_del/pow(@beta_min, \1.0)", line)
  line = re.sub(r"MB_MBMOD\((.*?)\)", r"MB_MODE && @mbin_mode == \1", line)
  line = re.sub(r"MB_BMAXCORR", r"MB_BMIN + MB_BDEL * MB_N", line)
  line = re.sub("MB_ISFIT", "@flags & MBEST_FIT", line)
  line = re.sub("SHK_GAUGE_VAL", "(@ens_invw[ia_] * @m) / (@it_et[ia_] - @is_et[ia_])", line)
  line = re.sub("SHK_FORMULA",  'amp t^(-exp)', line)
  line = re.sub("SHK_MODE_VALID", "@shk_mode >= 0 && @shk_mode <= 2", line)
  line = re.sub("SHK_MAX_VALID", "@shk_max < 0.9 && @shk_max >= 0.0", line)
  line = re.sub("SHK_MODE", "@shk_mode", line)
  line = re.sub("MB_BINIT", "0.5 * (MB_BMIN + MB_BMAX)", line)
  line = re.sub("MB_N_FORMULA", "(int)((MB_BMAX - MB_BMIN)/MB_BDEL - 1e-5) + 1", line)
  line = re.sub("MB_BMIN_VALID", "MB_BMIN >= 0.0", line)
  line = re.sub("MB_BMAX_VALID", "MB_BMAX > MB_BMIN", line)
  line = re.sub("MB_BARR_BAD", "@_check_beta_array( @ )", line)
  line = re.sub("MB_BDEL_VALID", "MB_BDEL > 1e-6", line)
  line = re.sub("MB_BETA_FUNC", "MB_BMIN + ia_ * MB_BDEL", line)
  line = re.sub("MB_BMAX", "@beta_max", line)
  line = re.sub("MB_BMIN", "@beta_min", line)
  line = re.sub("MB_BDEL", "@beta_del", line)
  line = re.sub("ENS_FORMULA", 'w = [1 + amp*exp^{-0.5 ((beta - focus)/dbeta)^2}] / beta^factor', line)
  line = re.sub("ENS_INVW_FUNC", "@_ens_invw_calc(@, 0.5*(@beta[ia_] + @beta[ia_+1]), NULL, NULL)", line)
  line = re.sub("ENS_INVWB_FUNC", "@_ens_invw_calc(@, @beta[ia_], NULL, NULL)", line)
  line = re.sub("ENS_FOCUS", "fabs(@ens_amp) > 1e-8", line)
  line = re.sub("MB_XSUMS_SIZE", "@m * @n", line)
  line = re.sub("MB_SBSIZE", "@sumsb_cnt * @n", line)
  line = re.sub("MB_HBSIZE", "@sumsb_cnt * (@n+1)", line)
  line = re.sub("MB_ORDER", "@order", line)
  line = re.sub("MB_MODE", "@mode", line)
  line = re.sub("MB_MCALC2", "@_m_calc_xsums( @ )", line)
  line = re.sub("MB_MCALC", "@_m_calc( @ )", line)
  line = re.sub("MB_SET_XSUMS", "@_setup_xsums(@, NULL)", line)
  line = re.sub("MB_MKWIN", "@_setup_windows( @ )", line)
  line = re.sub("MB_N", "@n", line)
  line = re.sub(r"MB_M([^\w]|$)", r"@m\1", line)
  return line
'''


r'''
file = "mb_t.def"
var_prefix = "eh_"
key_prefix = ""
flag_prefix = "MBEST_EHIST_"
xm_prefix = "EH_"
cond0_common = "@eh_mode"
def sub(line):
  line = re.sub(r"EH_MBDEL_VALID\((.*?)\)", r"@eh_mbin_delta > @beta_del/pow(MB_(beta_min), \1.0)", line)
  line = re.sub(r"EH_MBMOD\((.*?)\)", r"EH_MODE && @eh_mbin_mode == \1", line)
  line = re.sub('EH_CNT_FORMULA', '(int)((EH_MAX-EH_MIN)/EH_DEL-1e-5 + 1)', line)
  line = re.sub("EH_MIN", "@eh_min", line)
  line = re.sub("EH_MAX", "@eh_max", line)
  line = re.sub("EH_DEL", "@eh_del", line)
  line = re.sub("EH_NSTSAVE", "@eh_nstsave", line)
  line = re.sub("EH_SKIP_VALID", "@eh_skip > 0", line)
  line = re.sub("EH_MODE_VALID", "@eh_mode == 0 || @eh_mode == 1", line)
  line = re.sub("EH_MODE", "@eh_mode", line)
  line = re.sub("EHS_SIZE", "@eh_cnt * @n", line)
  line = re.sub("EHT_SIZE", "@eh_cnt", line)
  line = re.sub("MB_N", "@n", line)
  return line
'''


xm_item = xm_prefix + "ITEM"
xm_dumy = xm_prefix + "DUMY"

arg = r"([\s\w\.\+\@\-\*\/\=\>\<\!\[\]\{\}\(\)\&\|\"]+),"
sp_desc = r'\s*(\".*?(?<!\\)\"|\w+)'
pattern = r"\(" + arg * 9 + sp_desc + r"\s*\)"
print "pattern:", pattern


src = open(file).readlines()
block = []

cond0_once = 0
for line in src:
  if line.startswith(xm_item) or line.startswith(xm_dumy):
    m = re.search(pattern, line)
    if not m:
      print "cannot handle line:\n%s" % line
      exit(1)
    var   = m.group(1).strip()
    key   = m.group(2).strip()
    tp    = m.group(3).strip()
    sfx   = m.group(4).strip()
    cnt   = sub(m.group(5).strip())
    defl  = sub(m.group(6).strip())
    tag   = m.group(7).strip()
    cond0 = sub(m.group(8).strip())
    cond1 = sub(m.group(9).strip())
    desc  = sub(m.group(10).strip(' "'))
    cmds = ""
    
    if var != "" and tag != "FLAG": 
      var = var_prefix + var

    if tp.endswith("*"):
      tp = tp[:-1].strip()
      var = "*" + var

    if line.startswith(xm_dumy):
      if tag == "FLAG":
        cmds += "$flag: %-15s %s; " % (var, sfx)
      else:
        cmds += "$altvar: %s; " % (var)
      var = ""
      tp = ""
      
    if len(key):
      cmds += " $key: %s; " % (key_prefix + key)

    if cnt != "1":
      cmds += " $cnt: %s; " % cnt
    if defl != "" and not (tag == "FLAG" and key == ""):
      # remove ssdup
      pat1 = r'ssdup\(\"(.*?)\"\)' 
      m = re.match(pat1, defl);
      if m:
        defl = re.sub(pat1, r"\1", defl)
      cmds += " $def: %s; " % defl
    if tag == "MUST":
      if cond0 not in ("FALSE", "0", "NO"):
        cmds += " $key_must;"
      else:
        cmds += " $io:none;"
        cond0 = ""
    elif tag == "OPTN":
      pass
    elif tag == "FLAG":
      pass
    elif tag == "DARR" and cnt != "":
      pass
    elif tag == "SARR":
      var += sfx
    else:
      print "don't know how to handle tag = %s" % tag
      exit(1)

    if (cond0 != "" and cond0 not in ("1", "TRUE", "YES")):
      if cond0_common != "" and cond0 == cond0_common:
        if not cond0_once:
          cmds += " $prereq := %s; " % cond0
          cond0_once = 1
        else: pass
      else:
        cmds += " $prereq: %s; " % cond0
    if cond1 != "" and cond1 not in ("1", "TRUE", "YES"):
      cmds += " $valid: %s;" % cond1

    if len(desc) > 0 and desc[-1:] not in ";,." and cmds != "":
      desc += ";"
    if tag == "FLAG":
      comment = "/* " + cmds + " " + desc + " */"
    else:
      comment = "/* " + desc + cmds + " */"
    block += [(tp, var, comment)]
  else:
    if len(block) == 0: continue
    w1 = max(len(a[0]) for a in block)
    w2 = max(len(a[1]) for a in block) + 1
    for a in block:
      if a[0] != "":  # type var; comment
        print "  %-*s %-*s %s" % (w1, a[0], w2, a[1]+";", a[2])
      else:  # comment only
        print "  %s" % (a[2])
    print "  /* */"
    #raw_input()
    block = []
    
  
    

