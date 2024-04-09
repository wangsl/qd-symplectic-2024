/* $Id$ */

#ifndef LST_H
#define LST_H

#include <cassert>
#include <iostream>
using namespace std;
#include "die.h"
#include "indent.h"

template<class T> class LstItr;
template<class T> class LstPairItr;

template<class T> class Lst
{
  friend class LstItr<T>;
  friend class LstPairItr<T>;

protected:
  struct LstRep
  {
    T t;
    LstRep *next;
    LstRep(const T &s, LstRep *n = 0) : t(s), next(n) { }
  };
  LstRep *head, *tail;
  int n;
  
public:
  Lst() : head(0), tail(0), n(0) { }
  ~Lst() { remove_all(); }

  Lst<T> & operator=(const Lst<T> &lst)
  {
    remove_all();
    add(lst);
    return *this;
  }
  
  void remove_all()
  { 
    while (head) {
      tail = head->next;
      delete head;
      head = tail;
    }
    n = 0;
  }
  
  void remove_last()
  {
    for (LstRep *p = head; p && p->next; p = p->next)
      if (p->next == tail) {
	delete p->next;
	p->next = 0;
	tail = p;
	n--;
	return;
      }
  }
  
  void add(const T &t)
  { 
    n++;
    if (!head)
      tail = head = new LstRep(t);
    else
      tail = tail->next = new LstRep(t);
  }
  
  T & first() 
  { 
    assert(head);
    return head->t; 
  }
  
  const T & first() const 
  { 
    assert(head);
    return head->t; 
  }
  
  T & last() 
  {
    assert(tail);
    return tail->t;
  }
  
  const T & last() const
  {
    assert(tail);
    return tail->t;
  }
  
  int size() const { return n; }
  
  inline void add(const Lst<T> &l);
  inline Lst(const Lst<T> &);

  // add by Shenglong Wang
  inline void show_in_one_line(const char *header = 0) const;
};

template<class T> class LstItr
{
private:
  typename Lst<T>::LstRep *i;
public:
  LstItr() : i(0) { }
  void init(const Lst<T> &l) 
    { i = l.head; }
  int ok() const { return i != 0; }
  T & operator()() { return i->t; }
  void next() { i = i->next; }
};

template<class T> class LstPairItr
{
private:
  typename Lst<T>::LstRep *p, *q;
public:
  LstPairItr() : p(0), q(0) { }
  void init(const Lst<T> &l)
    { if ((p = l.head) != 0) q = p->next; }
  int ok() const { return p != 0 && q != 0; }
  T & i() { return p->t; }
  T & j() { return q->t; }
  void next() 
    { 
      if ((q = q->next) == 0 && (p = p->next) != 0) 
	q = p->next; 
    }
};

template<class T> inline void Lst<T>::add(const Lst<T> &l)
{
  LstItr<T> i;
  for (i.init(l); i.ok(); i.next())
    add(i());
}

template<class T> inline Lst<T>::Lst(const Lst<T> &l) : head(0), tail(0), n(0) 
{
  add(l);
}

template<class T> inline istream &operator>>(istream &s, Lst<T> &c)
{
  char bracket;
  s >> bracket;
  if (bracket != '{') { // Assume the list contains only one element
    s.putback(bracket);
    T t;
    s >> t;
    c.add(t);
    return s;
  }
  while(s) {
    s >> bracket;
    if (bracket == '}')
      return s;
    s.putback(bracket);

    if(bracket == '!' || bracket == '#') {
      char buf[1024];
      if(!s.getline(buf, 1024).good()) 
	die("Error reading class Lst, "
	    "line too long exceeded 1024 characters");
      continue;
    }
    
    T t;
    s >> t;
    c.add(t);
  }
  die("Lst: reached EOF before reading '}'");
  return s;
}

template<class T> inline ostream &operator<<(ostream &s, const Lst<T> &c)
{
  LstItr<T> i;
  s << "{\n";
  IndentPush();
  for (i.init(c); i.ok(); i.next())
    s << Indent() << i() << "\n";
  IndentPop();
  return s << Indent() << "}";
}

// add by Shenglong Wang

template<class T> inline void Lst<T>::show_in_one_line(const char *header) const
{
  if(header) cout << header;
  LstItr<T> i;
  for (i.init(*this); i.ok(); i.next())
    cout << " " << i();
  cout << endl;
}

#endif /* LST_H */
