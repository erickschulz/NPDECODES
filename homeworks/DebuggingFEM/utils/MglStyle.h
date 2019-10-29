#ifndef MGL_STYLE_H
#define MGL_STYLE_H

# include <deque>
# include <string>
# include <algorithm>
# include <mgl2/mgl.h>

/* create cross join B x A, such that:
 * A = {0, 1}, B = {a, b, bc}
 * => C = {0a, 1a, 0b, 1b, 0c, 1c}     */
static void crossjoin (const std::deque<std::string>& A, const std::deque<std::string>& B, std::deque<std::string>& res)
{
  for (auto b : B) {
    for (auto a : A) {
      res.push_back(a + b);
    }
  }
}

class MglStyle {
public:

  /* create new styles_ deque with styles from get_new */
  MglStyle() {
    get_new(styles_);
  }

  void get_new (std::deque<std::string>& new_deque) {
    std::deque<std::string> colors = { "b", "r", "g", "c", "m", "y", "G", "p", "o", "k", "n" },
      linetypes = { "-", ":", ";", "|", "j", "i", "=" };
      crossjoin(colors, linetypes, new_deque);
  }

  /* create new styles_ deque and remove 'already_used' */
  MglStyle (const std::string& already_used)
  {
    get_new(styles_);
    // if already_used is contained in styles_ remove it
    eliminate(already_used);
  }

  /* creates new styles_ deque and remove all strings in already_used */
  template <class Container>
  MglStyle (const Container& already_used_cont)
  {
    get_new(styles_);

    // iterate over all strings in already_used_cont and remove them
    for (auto already_used : already_used_cont) {
      eliminate(already_used);
    }

    // if all available styles have been used start from the beginning
    if (styles_.size() == 0) {
      get_new(styles_);
    }
  }

  std::string get_next ()
  {
    // if all available styles have been used start from the beginning
    if (styles_.size() == 0) {
      get_new(styles_);
    }
    std::string next = styles_[0];
    styles_.pop_front();
    return next;
  }

  void eliminate (const std::string& already_used) {
    auto it = std::find(styles_.begin(), styles_.end(), already_used);
    if (it != styles_.end()) {
      styles_.erase(it);
    }
  }

private:
  std::deque<std::string> styles_;
};

#endif
