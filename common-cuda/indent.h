

/* $Id$ */

#ifndef INDENT_H
#define INDENT_H

#ifdef __cplusplus
extern "C" {
#endif

  void IndentPush();
  void IndentPop();
  const char *Indent();

#ifdef __cplusplus
}
#endif

#endif /* INDENT_H */
