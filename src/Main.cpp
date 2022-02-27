/*
    Copyright (C) 2020
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is a part of APEX.

    APEX is distributed "AS IS" in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must
    be included in all copies or substantial portions of APEX.
*/
#include "Main.hpp"

std::string help_string =
  "\n"
  "  APEX: Toolkit for xQTL analysis\n"
  "     (c) 2019-2021 Corbin Quick and Li Guan.\n"
  "\n"
  "  Usage and options:\n"
  "     ./apex [mode] --help       Print help menu for [mode].\n"
  "\n"
  "  Preprocessing modes:\n"
  "     ./apex factor {OPTIONS}    Estimate latent factors from\n"
  "                                 molecular trait data.\n"
  "\n"
  "     ./apex lmm {OPTIONS}       Precompute terms for linear\n"
  "                                 mixed model analysis.\n"
  "\n"
  "  Analysis modes:\n"
  "     ./apex cis {OPTIONS}       Run cis-xQTL analysis.\n"
  "\n"
  "     ./apex trans {OPTIONS}     Run trans-xQTL analysis.\n"
  "\n"
  "  Meta-analysis modes:\n"
  "     ./apex store {OPTIONS}     Store vcov (LD) data for\n"
  "                                 xQTL meta-analysis or\n"
  "                                 data sharing.\n"
  "\n"
  "     ./apex meta {OPTIONS}      Single and multi-variant\n"
  "                                 xQTL meta-analysis from\n"
  "                                 sumstat and vcov files.\n"
  "\n"
  "  Contact: corbinq@gmail.com\n";

// ------------------------------------
//  Main function (determine running mode)
// ------------------------------------

int main(int argc, char* argv[])
{

// ------------------------------------
//  Hide cursor while running, restore on exit or interrupt
// ------------------------------------

  auto lam_kill =
    [] (int i) { restore_cursor(); std::cerr << "\nKilled.\n" << "\n"; exit(0); };

  signal(SIGINT, lam_kill);
  signal(SIGABRT, lam_kill);
  signal(SIGTERM, lam_kill);
  // signal(SIGTSTP, lam_wait);

#ifndef __APPLE__
  // This function is not implemented in OS X's stdlib.
	at_quick_exit (restore_cursor);
#endif
  atexit (restore_cursor);


// ------------------------------------
//  APEX main menu: Parse running mode.
// ------------------------------------

  std::unordered_map<std::string, mode_fun> map{
    {"cis", cis},
    {"trans", trans},
    {"factor", factor},
    {"lmm", lmm},
    {"meta", meta},
    {"store", store}
  };

  // Argument parsing using https://github.com/Taywee/args
  args::ArgumentParser p0("apex: GWAS/QTL Toolkit.", "Contact: corbinq@gmail.com.\n");
  args::Group apex_group(p0, "Options that affect all apex commands");
  args::Flag legacy_vcov(apex_group, "legacy_vcov", "Use legacy format for vcov files.", {"legacy-vcov"});
  args::HelpFlag help0(p0, "help", "Display this help menu", {'h', "help"});

  p0.Prog(argv[0]);

  std::string mode_summary = "\ncis: cis-xQTL analysis.\n\ntrans: trans-xQTL analysis.\nmeta: xQTL meta-analysis.\n\nstore: store xQTL vcov (LD) data.\n\n";

  args::MapPositional<std::string, mode_fun> mode(p0, "mode", mode_summary, map);
  mode.KickOut(true);

  const std::vector<std::string> args(argv + 1, argv + argc);

  try {
    auto next = p0.ParseArgs(args);
    bool use_legacy_vcov = args::get(legacy_vcov);
    global_opts::set_legacy_vcov(use_legacy_vcov);

    if (mode) {
      return args::get(mode)(argv[0], next, std::end(args));
    } else {
      std::cout << help_string;
    }
  }
  catch (args::Help) {
    std::cout << help_string;
    return 0;
  }
  catch (args::Error e) {
    std::cerr << "\nUnknown command line argument(s).\nPrinting help menu:\n" << std::endl;
    std::cerr << help_string;
    return 1;
  }
  return 0;
}