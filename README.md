# Dissimilarity Scoring of Bacterial Colony Patterns Using a Siamese Neural Network

Bacterial colony patterns exhibit a high degree of diversity, making it crucial to characterize them for a deeper understanding of the underlying biological processes. However, these patterns can be irregular, and traditional approaches used for more standardized branching patterns—such as those formed by mineral crystallization, snowflakes, etc.—are often inapplicable, despite visual similarities. Additionally, bacterial colony patterns can vary from batch to batch due to biological noise, although their qualitative features are generally preserved. Established metrics for photography (e.g. [Structural similarity index(SSIM)](https://en.wikipedia.org/wiki/Structural_similarity_index_measure)) often introduce significant errors and may not be reliable for analyzing these patterns.

# Siaemese Network
![Picture1](https://github.com/user-attachments/assets/6ebb1e9d-c2a8-4024-ab47-2a51cde33bf7)

To address this challenge, a Siamese network was developed to assess the similarities between different colony patterns. The Siamese network consists of two CNNs that compress colony images (either simulated or experimental) into smaller vectors, with a [contrastive loss](http://yann.lecun.com/exdb/publis/pdf/hadsell-chopra-lecun-06.pdf![image](https://github.com/user-attachments/assets/6e174647-4bd6-472c-b4bc-1695a9494b77)
) function applied to score their differences. During **training**, the CNN learns to extract key features from each class of patterns. The training data consists of simulated patterns from 5 distinct classes of branching patterns, characterized by varying combinations of branching width and density (in experiments, the patterns were generated under different environmental conditions; See this [paper](https://doi.org/10.15252/msb.202010089) for more details).

<p align="center">
Training Data
</p>


![Picture1 3 35 09 PM](https://github.com/user-attachments/assets/a894723f-4752-4a73-b522-ba853e6c3bcc)

In the **testing** phase, two images are fed into the trained model, which outputs a dissimilarity score. The results demonstrate that the model effectively scores different types of patterns: the more visually dissimilar the patterns are, the higher the dissimilarity score.

<p align="center">
Test Data
</p>
Both experimental (top row) and simulated (bottom row) data were used.


![Picture5](https://github.com/user-attachments/assets/2770b386-2ccc-4ee9-911c-24c12b2be92d)

# Examples
<p align="center">Dissimilarity score = 0.0.</p>
<p align="center">
  <img src="https://github.com/user-attachments/assets/92bcf2a4-87a0-4b67-922d-c54aeb601c85" alt="Picture2">
</p>

<p align="center">Dissimilarity score = 1.14.</p>
<p align="center">
  <img src="https://github.com/user-attachments/assets/e298be3a-1372-498d-985e-b28f5724cdaa" alt="Picture3">
</p>

<p align="center">Dissimilarity score = 1.53.</p>
<p align="center">
  <img src="https://github.com/user-attachments/assets/568399bd-ae4d-413e-95ed-b984b7caf092" alt="Picture4">
</p>
