#!/usr/bin/env ruby

# Input: clusters printed by Hammer-IT
# Output: FASTA file with kmers from these clusters
#         Name format: [cluster index].[kmer index in the cluster]-cov_[kmer count]-qual_[kmer quality]

FMT = ">%06d.%03d-cov_%d-qual_%f"

KMer = Struct.new :sequence, :count, :quality
Cluster = Struct.new :index, :kmers

clusters = []

data = File.read 'clusters.txt'
data.scan(/(\d+): \{ \n(.*?)center:.*?\}/m).each do |m|
  kmers = []
  m[1].split(", \n").each do |k|
    kmer_data = k.match(/(\w+): \((\d+), ([.0-9]+)\)/).to_a
    kmers << KMer.new(kmer_data[1], kmer_data[2].to_i, kmer_data[3].to_f)
  end
  clusters << Cluster.new(m[0].to_i, kmers)
end

clusters.each do |cluster|
  modes = cluster.kmers.sort_by(&:count).reverse
  modes.each_with_index do |mode, i|
    puts (FMT % [cluster.index, i, mode.count, mode.quality])
    puts mode.sequence
  end
end
