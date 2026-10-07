#ifndef TOP_CROSS_SECTION_LUMI_MASK_H
#define TOP_CROSS_SECTION_LUMI_MASK_H

#include <fstream>
#include <map>
#include <regex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class LumiMask {
 public:
  explicit LumiMask(const std::string &jsonPath) { Load(jsonPath); }

  bool Accept(unsigned int run, unsigned int luminosityBlock) const {
    const auto runEntry = intervals_.find(run);
    if (runEntry == intervals_.end()) return false;
    for (const auto &interval : runEntry->second)
      if (luminosityBlock >= interval.first &&
          luminosityBlock <= interval.second)
        return true;
    return false;
  }

  std::size_t NumberOfRuns() const { return intervals_.size(); }

 private:
  void Load(const std::string &jsonPath) {
    std::ifstream input(jsonPath);
    if (!input) throw std::runtime_error("Cannot open lumi mask: " + jsonPath);
    const std::string content((std::istreambuf_iterator<char>(input)),
                              std::istreambuf_iterator<char>());

    std::size_t position = 0;
    const std::regex pairPattern(
        R"(\[\s*([0-9]+)\s*,\s*([0-9]+)\s*\])");
    while ((position = content.find('"', position)) != std::string::npos) {
      const std::size_t quoteEnd = content.find('"', position + 1);
      if (quoteEnd == std::string::npos) break;
      const std::string runText =
          content.substr(position + 1, quoteEnd - position - 1);
      if (runText.empty() ||
          runText.find_first_not_of("0123456789") != std::string::npos) {
        position = quoteEnd + 1;
        continue;
      }

      const std::size_t listStart = content.find('[', quoteEnd);
      if (listStart == std::string::npos) break;
      int depth = 0;
      std::size_t listEnd = listStart;
      for (; listEnd < content.size(); ++listEnd) {
        if (content[listEnd] == '[') ++depth;
        if (content[listEnd] == ']' && --depth == 0) break;
      }
      if (listEnd >= content.size())
        throw std::runtime_error("Malformed lumi mask JSON");

      const unsigned int run = static_cast<unsigned int>(std::stoul(runText));
      const std::string list =
          content.substr(listStart, listEnd - listStart + 1);
      for (std::sregex_iterator match(list.begin(), list.end(), pairPattern), end;
           match != end; ++match) {
        intervals_[run].emplace_back(
            static_cast<unsigned int>(std::stoul((*match)[1].str())),
            static_cast<unsigned int>(std::stoul((*match)[2].str())));
      }
      position = listEnd + 1;
    }
    if (intervals_.empty())
      throw std::runtime_error("No run intervals found in lumi mask: " + jsonPath);
  }

  std::map<unsigned int,
           std::vector<std::pair<unsigned int, unsigned int>>> intervals_;
};

#endif

