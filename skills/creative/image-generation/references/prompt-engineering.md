# Prompt Engineering for Image Generation

## Overview

Effective prompt engineering is the key to generating high-quality images with AI models. This guide covers proven techniques, patterns, and best practices for crafting prompts that produce consistent, high-quality results across both Gemini and DALL-E models.

## Fundamental Principles

### 1. Be Specific

Vague prompts produce unpredictable results. Specificity guides the model toward your vision.

**Bad:**
```
A nice landscape
```

**Good:**
```
Mountain valley at sunrise, mist rolling over alpine lake, pine trees in
foreground, warm golden light, peaceful atmosphere, photorealistic landscape
photography
```

**Why it's better:** Includes subject (mountain valley), time (sunrise), elements (lake, trees, mist), lighting (warm golden), mood (peaceful), and style (photorealistic).

### 2. Structure Your Prompts

Follow a consistent structure to ensure all important elements are included:

**Basic Structure:**
```
[Subject] + [Style] + [Details] + [Lighting] + [Mood/Atmosphere] + [Technical Details]
```

**Example:**
```
Subject: Portrait of a jazz musician
Style: Film noir photography
Details: Playing saxophone, wearing fedora hat, smoke in air
Lighting: Dramatic side lighting, high contrast
Mood: Moody and atmospheric
Technical: Black and white, shot on 50mm lens, shallow depth of field
```

**Full Prompt:**
```
Portrait of a jazz musician playing saxophone, wearing fedora hat with smoke in
the air, film noir photography style, dramatic side lighting with high contrast,
moody and atmospheric, black and white, shot on 50mm lens, shallow depth of field
```

### 3. Specify Style Clearly

Style descriptors dramatically affect the output. Be explicit about the artistic style you want.

**Photography Styles:**
- `photorealistic`, `professional photography`, `commercial photography`
- `portrait photography`, `landscape photography`, `macro photography`
- `street photography`, `documentary style`, `editorial photography`
- `film noir`, `vintage photography`, `instant camera polaroid`

**Artistic Styles:**
- `digital art`, `digital painting`, `concept art`
- `oil painting`, `watercolor`, `acrylic painting`
- `pencil sketch`, `charcoal drawing`, `ink illustration`
- `vector art`, `flat design`, `isometric illustration`

**3D and Rendering:**
- `3D render`, `octane render`, `unreal engine`
- `photorealistic 3D`, `stylized 3D`, `low poly 3D`

**Animation and Comics:**
- `anime style`, `manga style`, `cartoon illustration`
- `Disney style`, `Studio Ghibli style`, `pixel art`
- `comic book art`, `graphic novel style`

### 4. Include Lighting Information

Lighting transforms the mood and quality of images.

**Natural Lighting:**
- `golden hour lighting`, `blue hour`, `midday sun`
- `overcast soft light`, `sunrise`, `sunset`
- `dappled sunlight through trees`, `window light`

**Studio Lighting:**
- `studio lighting`, `three-point lighting`, `softbox lighting`
- `rim lighting`, `backlit`, `silhouette lighting`
- `key light from left`, `fill light`, `hair light`

**Dramatic Lighting:**
- `dramatic lighting`, `chiaroscuro`, `low key lighting`
- `high contrast`, `side lighting`, `dramatic shadows`
- `spotlight`, `stage lighting`, `moody lighting`

**Ambient and Special:**
- `soft diffused light`, `ambient lighting`, `natural light`
- `neon lighting`, `bioluminescent glow`, `candlelight`
- `volumetric lighting`, `god rays`, `atmospheric lighting`

### 5. Specify Composition

Guide how elements are arranged in the frame.

**Camera Angles:**
- `eye level`, `low angle looking up`, `high angle looking down`
- `bird's eye view`, `overhead shot`, `worm's eye view`
- `dutch angle`, `over the shoulder`, `point of view shot`

**Framing:**
- `close-up`, `medium shot`, `wide shot`, `extreme close-up`
- `full body shot`, `portrait`, `headshot`
- `centered composition`, `rule of thirds`, `symmetrical composition`

**Depth:**
- `shallow depth of field`, `bokeh background`
- `foreground interest with blurred background`
- `everything in sharp focus`, `deep focus`

## Advanced Techniques

### Layered Descriptions

Build complexity by describing foreground, midground, and background separately.

```
Foreground: A vintage typewriter with paper inserted, hands typing
Midground: Scattered manuscript pages and coffee cup on wooden desk
Background: Bookshelf filled with classic literature, soft window light
Style: Warm cozy aesthetic, film photography, nostalgic mood
```

### Negative Space

Describe what should NOT be in the image to guide composition.

```
Minimalist product photography of wireless earbuds, vast white negative space
surrounding subject, clean and spacious composition, professional studio shot
```

### Material and Texture Specification

Detail the physical properties of objects.

**Materials:**
- `brushed metal`, `polished chrome`, `matte black finish`
- `soft fabric`, `rough weathered wood`, `smooth glass`
- `translucent plastic`, `glossy ceramic`, `textured leather`

**Textures:**
- `fine grain texture`, `rough surface`, `smooth polished`
- `weathered and aged`, `pristine and new`, `worn and vintage`

### Emotional and Atmospheric Descriptors

Add mood and feeling to guide the overall atmosphere.

**Moods:**
- `serene and peaceful`, `energetic and dynamic`, `melancholic and contemplative`
- `joyful and vibrant`, `mysterious and enigmatic`, `dramatic and intense`

**Atmospheres:**
- `cozy and intimate`, `vast and expansive`, `claustrophobic and tense`
- `dreamlike and surreal`, `realistic and grounded`, `magical and whimsical`

## Model-Specific Strategies

### Gemini Models

**Gemini 2.5 Flash:**
- Keep prompts focused but detailed
- Good for iterative refinement via chat
- Responds well to clear style directions

**Gemini 3 Pro:**
- Can handle very complex, multi-element descriptions
- Better at following detailed instructions
- Excellent for text-in-image (specify exact text)
- Use multiple reference images when needed

**Example Gemini Prompt:**
```
Create a professional tech conference poster featuring "AI Summit 2025" as the
main title in bold modern sans-serif font. Background: abstract neural network
visualization with glowing blue and purple nodes and connections. Foreground:
silhouettes of diverse professionals networking. Layout: title at top third,
date "March 15-17" in smaller text, clean modern design, professional gradient
from deep blue to purple, high resolution poster format
```

### DALL-E Models

**DALL-E 3:**
- Use natural, conversational language
- The model expands prompts automatically
- Focus on essence and key elements
- Trust the model's interpretation
- Good for creative, artistic requests

**Example DALL-E 3 Prompt:**
```
A cozy independent bookstore interior with comfortable reading nooks, warm
lighting from vintage lamps, customers browsing shelves and reading in armchairs,
large window showing rainy street outside, plants on windowsills, artistic and
inviting atmosphere
```

**DALL-E 2:**
- Be very explicit and literal
- Include all details you want
- Less interpretation, more direct
- Good for controlled, predictable outputs

**Example DALL-E 2 Prompt:**
```
Interior of independent bookstore, wooden bookshelves along walls filled with
books, two armchairs with reading lamps, three customers browsing, large window
on left wall showing rainy street, potted plants on windowsill, warm yellow
lighting, artistic digital art style, cozy atmosphere
```

## Common Use Case Patterns

### Product Photography

```
[Product name/description], professional product photography, white seamless
background, studio lighting with soft shadows, [angle/orientation], clean
minimal composition, high resolution, sharp focus, commercial quality,
[specific color/material details]
```

**Example:**
```
Premium leather wallet in cognac brown, professional product photography, white
seamless background, studio lighting with soft shadows, 3/4 angle showing
interior card slots, clean minimal composition, high resolution, sharp focus,
commercial quality, visible grain texture
```

### Character Design

```
[Character description] character design, [body type/pose], [clothing/accessories],
[artistic style], [color palette], [mood/personality], [background context],
detailed digital art, character sheet format
```

**Example:**
```
Female cyberpunk hacker character design, athletic build in confident standing
pose, tactical jacket with LED strips and cargo pants, techwear style, neon blue
and black color palette with orange accents, confident and rebellious mood,
subtle futuristic city background, detailed digital art, character concept art
```

### Logo Design

```
[Company/concept] logo, [style] design, [geometric/symbolic elements],
[color scheme], [mood/feeling], professional, [format specifications],
suitable for [use cases]
```

**Example:**
```
Sustainable energy company logo, minimalist geometric design, stylized sun with
leaf elements forming circular shape, green and gold color scheme, modern and
trustworthy feeling, professional, vector style, clean lines, suitable for
digital and print applications
```

### Landscape/Environment

```
[Location/scene type], [time of day], [weather/atmospheric conditions],
[key natural elements], [foreground/midground/background], [mood],
[photography/art style], [technical specs]
```

**Example:**
```
Norwegian fjord landscape, late evening during blue hour, light fog over water,
jagged mountain peaks in background, lone wooden cabin on shoreline in
foreground, serene and mystical mood, landscape photography style, shot on
24mm wide angle lens, high detail
```

### Food Photography

```
[Dish description], professional food photography, [plating/presentation],
[background/context], [lighting], [angle], appetizing and fresh, high quality
restaurant photography, [garnishes/styling]
```

**Example:**
```
Gourmet ramen bowl, professional food photography, artfully plated with soft-
boiled egg halves, pork chashu, green onions, and nori, rustic wooden table
background with chopsticks, natural window lighting with soft shadows, 45-degree
overhead angle, appetizing and fresh, high quality restaurant photography,
steam rising from hot broth
```

### Architectural Visualization

```
[Building/space type], [architectural style], [exterior/interior], [time of day],
[materials and textures], [context/surroundings], [mood], architectural
photography/rendering, [technical specs]
```

**Example:**
```
Modern minimalist residential home exterior, contemporary architecture with
clean lines and large glass panels, during golden hour sunset, white stucco
walls with natural wood accents, landscaped garden with native plants,
peaceful suburban setting, architectural photography, shot with tilt-shift
lens, high detail
```

## Iteration Strategies

### Progressive Refinement

Start broad, then add detail:

**Iteration 1:**
```
A forest scene with a cabin
```

**Iteration 2:**
```
A log cabin in a pine forest, surrounded by tall trees, forest clearing
```

**Iteration 3:**
```
A rustic log cabin in a pine forest clearing, surrounded by tall evergreen
trees, sunset lighting filtering through branches, smoke from chimney, cozy
atmosphere, landscape photography style
```

### Style Exploration

Keep subject constant, vary style:

**Base Subject:** "A cat sitting on a windowsill"

**Photorealistic:**
```
A cat sitting on a windowsill, professional pet photography, natural window
light, soft focus background, high detail fur texture
```

**Illustration:**
```
A cat sitting on a windowsill, children's book illustration style, watercolor
technique, soft pastel colors, gentle and whimsical
```

**Artistic:**
```
A cat sitting on a windowsill, impressionist oil painting style, visible
brushstrokes, warm color palette, artistic interpretation
```

### Mood Variations

Keep subject and style, vary atmosphere:

**Base:** "City street at night, cinematic photography"

**Moody:**
```
City street at night, cinematic photography, rain-slicked pavement reflecting
neon signs, moody and atmospheric, film noir aesthetic
```

**Vibrant:**
```
City street at night, cinematic photography, vibrant neon lights, bustling
energy, colorful and lively atmosphere
```

**Desolate:**
```
City street at night, cinematic photography, empty and quiet, single streetlight,
mysterious and isolated mood, urban solitude
```

## Common Mistakes and Fixes

### Mistake 1: Too Vague

**Problem:**
```
A beautiful sunset
```

**Fix:**
```
Coastal sunset over rocky beach, dramatic orange and purple clouds, waves
crashing on shore, silhouetted palm tree in foreground, warm golden hour
lighting, landscape photography
```

### Mistake 2: Conflicting Styles

**Problem:**
```
Photorealistic cartoon character in anime style oil painting
```

**Fix (choose one style):**
```
Anime character portrait, digital anime art style, detailed shading, vibrant
colors, professional anime illustration
```

### Mistake 3: Missing Technical Details

**Problem:**
```
Portrait of a person
```

**Fix:**
```
Portrait of a middle-aged woman, professional portrait photography, natural
expression, soft studio lighting, neutral gray background, shot on 85mm lens,
shallow depth of field, high quality
```

### Mistake 4: Overcomplication

**Problem:**
```
A hyperrealistic 8K ultra-detailed photorealistic quantum-rendered supremely
detailed hyper-realistic ultra-HD masterpiece...
```

**Fix:**
```
Portrait photography, professional quality, high detail, 85mm lens
```

### Mistake 5: Unclear Composition

**Problem:**
```
A cat, a dog, a bird, and a fish together
```

**Fix:**
```
Indoor pet gathering: cat and dog sitting together in foreground, birdcage on
table in midground, fish tank visible in background, cozy living room setting,
natural window light
```

## Testing and Validation

### A/B Testing Prompts

Generate with slight variations to understand what works:

**Version A:**
```
Futuristic city skyline, cyberpunk aesthetic, neon lights
```

**Version B:**
```
Futuristic city skyline at night, cyberpunk aesthetic with prominent neon
signs and holographic advertisements, rain-slicked streets, atmospheric fog
```

Compare results to learn which elements improve output quality.

### Prompt Journaling

Keep a record of successful prompts:

```markdown
## Successful Prompts Log

### Character Portraits
- **Style:** "detailed digital art, painterly style, soft brush strokes"
- **Lighting:** "soft side lighting with gentle rim light"
- **Quality:** "high detail, professional illustration"

### Product Photos
- **Background:** "white seamless background" OR "natural lifestyle setting"
- **Lighting:** "studio lighting with soft shadows"
- **Angle:** "3/4 angle view" works best for most products
```

## Quick Reference Formulas

### Photorealistic Images
```
[Subject] + professional photography + [lighting type] + [camera/lens] +
[composition] + high detail
```

### Illustrations
```
[Subject] + [art style] illustration + [color palette] + [mood] +
[technical style markers]
```

### Logos
```
[Concept] logo + [style] design + [shapes/symbols] + [colors] +
professional + vector style
```

### Scenes/Environments
```
[Location] + [time/weather] + [key elements] + [mood] + [style] +
[camera specs]
```

## Resources and Inspiration

- Study successful AI art on platforms like Midjourney showcase, DALL-E gallery
- Analyze prompts that produced great results
- Learn photography terminology to describe lighting/composition
- Build a personal style guide with proven prompt patterns
- Experiment regularly and document results

## Summary Checklist

Before submitting a prompt, ensure you've included:

- [ ] Clear subject/main focus
- [ ] Specific style (photography, illustration, 3D, etc.)
- [ ] Lighting description
- [ ] Mood/atmosphere
- [ ] Composition/camera angle (if applicable)
- [ ] Important details (colors, textures, elements)
- [ ] Technical specifications (lens, quality markers)
- [ ] Context/setting/background

The more intentional and specific your prompts, the more consistent and high-quality your generated images will be!
