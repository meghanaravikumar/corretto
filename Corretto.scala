/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.bdgenomics.corretto

import org.apache.spark.SparkContext._
import org.apache.spark.{ Logging, SparkContext }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.corretto.rdd.ErrorCorrection
import org.bdgenomics.formats.avro._
import org.bdgenomics.utils.cli._
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }

object Corretto extends BDGCommandCompanion {
  val commandName = "corretto"
  val commandDescription = "Correct errors in reads"

  def apply(cmdLine: Array[String]) = {
    new Corretto(Args4j[CorrettoArgs](cmdLine))
  }
}

class CorrettoArgs extends Args4jBase with ParquetArgs {
  @Argument(required = true, metaVar = "INPUT", usage = "The reads to correct.", index = 0)
  var reads: String = null

  @Argument(required = true, metaVar = "OUTPUT", usage = "The corrected reads.", index = 1)
  var outputPath: String = null

  @Args4jOption(required = false, name = "-kmer_length", usage = "The k-mer length to use.")
  var kmerLen: Int = 20

  @Args4jOption(required = false, name = "-maxIterations", usage = "The maximum number of iterations to run.")
  var maxIterations: Int = 20

  @Args4jOption(required = false, name = "-fixingThreshold", usage = "The quality score threshold to fix below.")
  var threshold: Int = 20

  @Args4jOption(required = false, name = "-errorRate", usage = "The baseline expected error rate.")
  var errorRate: Double = 0.05

  @Args4jOption(required = false, name = "-ploidy", usage = "The average sample ploidy.")
  var ploidy: Int = 2
}

class Corretto(protected val args: CorrettoArgs) extends BDGSparkCommand[CorrettoArgs] {
  val companion = Corretto

  def run(sc: SparkContext) {
    // load reads
    val reads = sc.loadAlignments(args.reads)

    // correct reads
    val correctedReads = ErrorCorrection(reads,
                                         args.kmerLen,
                                         args.maxIterations,
                                         args.threshold,
                                         args.errorRate,
                                         args.ploidy)
    // save
    correctedReads.adamParquetSave(args.outputPath)
  }
}

