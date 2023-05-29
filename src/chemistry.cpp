// Copyright 2020, the Aether Development Team (see doc/dev_team.md for members)
// Full license can be found in License.md

#include <string>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <random>

#include "aether.h"

using hashtable = std::unordered_map<std::string, int>;

// -----------------------------------------------------------------------------
// Initialize chemistry class
// -----------------------------------------------------------------------------

Chemistry::Chemistry(Neutrals neutrals,
                     Ions ions,
                     Inputs args,
                     Report &report) {

  std::string function = "Chemistry::Chemistry";
  static int iFunction = -1;
  report.enter(function, iFunction);

  read_chemistry_file(neutrals, ions, args, report);

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// Read chemistry CSV file
// -----------------------------------------------------------------------------

int Chemistry::read_chemistry_file(Neutrals neutrals,
                                   Ions ions,
                                   Inputs args,
                                   Report &report) {

  std::string function = "Chemistry::read_chemistry_file";
  static int iFunction = -1;
  report.enter(function, iFunction);

  std::ifstream infile_ptr;
  int iErr = 0;
  reaction_type reaction;

  report.print(1, "Reading Chemistry File : " + args.get_chemistry_file());

  infile_ptr.open(args.get_chemistry_file());

  if (!infile_ptr.is_open()) {
    std::cout << "Could not open chemistry file!\n";
    iErr = 1;
  } else {

    if (infile_ptr.good()) {

      std::vector<std::vector<std::string>> csv = read_csv(infile_ptr);

      int nLines = csv.size();

      if (nLines <= 2)
        iErr = 1;

      else {

        nReactions = 0;
        
        // Hash Map storing headers
        hashtable headers;
        for(int i = 0; i < csv[0].size(); ++i)
            headers[csv[0][i]] = i;

        // Skip 2 lines of headers!
        for (int iLine = 2; iLine < nLines; iLine++) {
          // Some final rows can have comments in them, so we want to
          // skip anything where the length of the string in column 2
          // is == 0:
          if (csv[iLine][headers["rate"]].length() > 0) {
            report.print(3, "interpreting chemistry line : " + csv[iLine][headers["loss1"]]);
            reaction = interpret_reaction_line(neutrals, ions,
                                               csv[iLine], report);
          }

          // check if it is part of a piecewise function,
          //   if so use sources/losses for last reaction
          if (reaction.nLosses == 0 && reaction.nSources == 0) {
            reaction.sources_names = reactions.back().sources_names;
            reaction.losses_names = reactions.back().losses_names;

            reaction.sources_ids = reactions.back().sources_ids;
            reaction.losses_ids = reactions.back().losses_ids;

            reaction.sources_IsNeutral = reactions.back().sources_IsNeutral;
            reaction.losses_IsNeutral = reactions.back().losses_IsNeutral;

            reaction.nLosses = reactions.back().nLosses;
            reaction.nSources = reactions.back().nSources;

            reaction.branching_ratio = reactions.back().branching_ratio;

            reaction.energy = reactions.back().energy;

            reaction.piecewiseVar = reactions.back().piecewiseVar;
          }

          if (reaction.nLosses > 0 && reaction.nSources > 0) {
            if (report.test_verbose(3))
              display_reaction(reaction);

            reactions.push_back(reaction);
            nReactions++;
          }
        }
          
        // add perturb to the all the reaction lines
        if(headers.find("perturb") != headers.end()) {
          json values = args.get_perturb_values();
            if(values.size() > 0) {
              if(values[0] == "all") {
                for(int i = 0; i < reactions.size(); ++i) {
                  precision_t perturb_rate = stoi(csv[i+2][headers["perturb"]]);
                  precision_t stdv = 0;
                  if(perturb_rate.length() > 0)
                      precision_t stdv = perturb_rate * reactions[i].rate;
                  std::random_device rd {};
                  std::mt19937 gen {rd()};
                  std::normal_distribution<> d {reactions[i], stdv};
                  reactions[i].rate = d(gen);
                }
              }
              else {
                for(auto &i : values) {
                  int line = stoi(i.substr(1));
                  precision_t perturb_rate = stoi(csv[i+2][headers["perturb"]]);
                  precision_t stdv = 0;
                  if(perturb_rate.length() > 0)
                    precision_t stdv = perturb_rate * reactions[i].rate;
                  std::random_device rd {};
                  std::mt19937 gen {rd()};
                  std::normal_distribution<> d {reactions[i], stdv};
                  reactions[i].rate = d(gen);
                }
              }
            }
          }
        }
      }
    }
  }

  report.exit(function);
  return iErr;
}

// -----------------------------------------------------------------------------
// Interpret a comma separated line of the chemical reaction file
// -----------------------------------------------------------------------------

Chemistry::reaction_type Chemistry::interpret_reaction_line(Neutrals neutrals,
                                                            Ions ions,
                                                            std::vector<std::string> line,
                                                            Report &report,
                                                            hashtable &headers) {

  std::string function = "Chemistry::interpret_reaction_line";
  static int iFunction = -1;
  report.enter(function, iFunction);

  reaction_type reaction;

  int i;
  int id_;
  bool IsNeutral;

  // Losses (left side) first:
  reaction.nLosses = 0;

  for (i = headers["loss1"]; i <= headers["loss3"]; i++) {
    find_species_id(line[i], neutrals, ions, id_, IsNeutral, report);

    if (id_ >= 0) {
      reaction.losses_names.push_back(line[i]);
      reaction.losses_ids.push_back(id_);
      reaction.losses_IsNeutral.push_back(IsNeutral);
      reaction.nLosses++;
    }
  }

  // Sources (right side) second:
  reaction.nSources = 0;

  for (i = headers["source1"]; i <= headers["source3"]; i++) {
    find_species_id(line[i], neutrals, ions, id_, IsNeutral, report);

    if (id_ >= 0) {
      reaction.sources_names.push_back(line[i]);
      reaction.sources_ids.push_back(id_);
      reaction.sources_IsNeutral.push_back(IsNeutral);
      reaction.nSources++;
    }
  }

  // Reaction Rate:
  reaction.rate = stof(line[headers["rate"]]);

  // Branching Ratio:
  if (line[headers["branching"]].length() > 0)
    reaction.branching_ratio = stof(line[headers["branching"]]);

  else
    reaction.branching_ratio = 1;


  // energy released as exo-thermic reaction:
  if (line[headers["heat"]].length() > 0)
    reaction.energy = stof(line[headers["heat"]]);

  else
    reaction.energy = 0;

  // default to zero (no piecewise, no exponent)
  reaction.min = 0;
  reaction.max = 0;
  reaction.type = 0;

  // if richards, check for temperature dependence
  if (headers["branching"] = 10) {
    //std::cout << line[17] << ", " << line[18] << ", " << line[19] << "\n";
    if (line[headers["Numerator"]].length() > 0) {
      reaction.numerator = stof(line[headers["Numerator"]]);
      reaction.denominator =    line[headers["Denominator"]];

      if (line[headers["Exponent"]].length() > 0)
        reaction.exponent = stof(line[headers["Exponent"]]);
    } else {
      // default to 0 (calc_chemical_sources will use constant rate)
      reaction.type = 0;
    }

    reaction.piecewiseVar =      line[headers["Piecewise"]];

    //std::cout << line[10] << ", " << line[17] << ", " << line[17+4] << "\n";
    if (line[headers["Min"]].length() > 0)
      reaction.min =        stoi(line[headers["Min"]);

    if (line[headers["Max"].length() > 0)
      reaction.max =        stoi(line[headers["Max"]);

    if (line[headers["Formula Type"].length() > 0)
      reaction.type =       stoi(line[headers["Formula Type"]]);
  }

  report.exit(function);
  return reaction;
}

// -----------------------------------------------------------------------------
// Match a string to the neutral or ion species
// -----------------------------------------------------------------------------

void Chemistry::find_species_id(std::string name,
                                Neutrals neutrals,
                                Ions ions,
                                int &id_,
                                bool &IsNeutral,
                                Report &report) {

  std::string function = "Chemistry::find_species_id";
  static int iFunction = -1;
  report.enter(function, iFunction);

  int iSpecies;
  IsNeutral = false;

  id_ = neutrals.get_species_id(name, report);

  if (id_ > -1)
    IsNeutral = true;

  else
    id_ = ions.get_species_id(name, report);

  report.exit(function);
  return;
}

// -----------------------------------------------------------------------------
// Display a reaction:
// -----------------------------------------------------------------------------

void Chemistry::display_reaction(Chemistry::reaction_type reaction) {

  int i;

  std::cout << "Number of Losses : " << reaction.nLosses << "\n";
  std::cout << "Number of Sources : " << reaction.nSources << "\n";

  for (i = 0; i < reaction.nLosses; i++)
    std::cout << reaction.losses_names[i] << " + ";

  std::cout << " -> ";

  for (i = 0; i < reaction.nSources; i++)
    std::cout << reaction.sources_names[i] << " + ";

  std::cout << " ( RR : " << reaction.rate << ")\n";

  for (i = 0; i < reaction.nLosses; i++)
    std::cout << reaction.losses_ids[i]
              << "(" << reaction.losses_IsNeutral[i] << ")" << " + ";

  std::cout << " -> ";

  for (i = 0; i < reaction.nSources; i++)
    std::cout << reaction.sources_ids[i]
              << "(" << reaction.sources_IsNeutral[i]
              << ")" << " + ";

  std::cout << " ( RR : " << reaction.rate << ")\n";

  if (reaction.type > 0) {
    std::cout << "Temperature Dependence: ("
              << reaction.numerator
              << "/"
              << reaction.denominator
              << ")^"
              << reaction.exponent << "\n";


  }

  if (reaction.min < reaction.max) {
    std::cout << "Range: "
              << reaction.min
              << " < "
              << reaction.piecewiseVar;

    if (reaction.max)
      std::cout << " < " << reaction.max;

    std::cout << "\n";
  }
}
