# Text-to-Image Generation Examples

> **Note**: Examples use the `google-genai` SDK. The `google-generativeai` package has been deprecated.

## Overview

Text-to-image generation is the core capability of AI image models. This guide provides practical examples across different use cases, showing effective prompts and expected workflows.

## Basic Workflow

1. Understand user requirements
2. Select appropriate model
3. Craft detailed prompt
4. Generate image
5. Save and present result
6. Offer refinements if needed

## Example 1: Photorealistic Portrait

### User Request
"Generate a professional headshot for a LinkedIn profile"

### Analysis
- **Type:** Portrait photography
- **Style:** Professional, photorealistic
- **Use case:** LinkedIn (square format ideal)
- **Model choice:** Gemini 2.5 Flash (fast, adequate quality) or DALL-E 3 (higher quality)

### Prompt
```
Professional business headshot portrait, middle-aged professional in business
casual attire, confident friendly expression, neutral gray background, soft
studio lighting, shot on 85mm lens, shallow depth of field, high quality
professional photography
```

### Code Implementation (Gemini)
```python
from google import genai
from google.genai import types

client = genai.Client()  # Uses GEMINI_API_KEY env var

prompt = """Professional business headshot portrait, middle-aged professional
in business casual attire, confident friendly expression, neutral gray
background, soft studio lighting, shot on 85mm lens, shallow depth of field,
high quality professional photography"""

response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=prompt,
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        image_config=types.ImageConfig(aspect_ratio="1:1")
    )
)

for part in response.parts:
    if part.inline_data is not None:
        image = part.as_image()
        image.save("linkedin_headshot.png")
        print("Generated: linkedin_headshot.png")
```

### Code Implementation (DALL-E 3)
```python
from openai import OpenAI
from pathlib import Path
import requests

client = OpenAI(api_key=os.environ["OPENAI_API_KEY"])

prompt = """Professional business headshot portrait, confident professional
in business casual attire, friendly expression, neutral gray background,
soft studio lighting, high quality photography"""

response = client.images.generate(
    model="dall-e-3",
    prompt=prompt,
    size="1024x1024",
    quality="hd"
)

image_url = response.data[0].url
image_data = requests.get(image_url).content
Path("linkedin_headshot.png").write_bytes(image_data)
print("Generated: linkedin_headshot.png")
```

## Example 2: Landscape Photography

### User Request
"Create a scenic mountain landscape for a travel blog header"

### Analysis
- **Type:** Landscape photography
- **Style:** Photorealistic, scenic
- **Use case:** Blog header (wide format)
- **Model choice:** Gemini 2.5 Flash or DALL-E 3 (landscape orientation)

### Prompt
```
Majestic mountain landscape at golden hour, snow-capped peaks under warm sunset
light, alpine valley with turquoise glacial lake in foreground, scattered pine
trees, wispy clouds, serene and breathtaking atmosphere, professional landscape
photography, shot on 24mm wide angle lens, high detail
```

### Code Implementation (Gemini)
```python
prompt = """Majestic mountain landscape at golden hour, snow-capped peaks under
warm sunset light, alpine valley with turquoise glacial lake in foreground,
scattered pine trees, wispy clouds, serene and breathtaking atmosphere,
professional landscape photography, shot on 24mm wide angle lens, high detail"""

response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents=prompt,
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        image_config=types.ImageConfig(aspect_ratio="16:9")  # Wide format for header
    )
)

for part in response.parts:
    if part.inline_data is not None:
        image = part.as_image()
        image.save("travel_blog_header.png")
```

### Code Implementation (DALL-E 3)
```python
response = client.images.generate(
    model="dall-e-3",
    prompt="""Majestic mountain landscape at golden hour, snow-capped peaks,
    alpine valley with turquoise lake, pine trees, breathtaking atmosphere,
    landscape photography""",
    size="1792x1024",  # Landscape orientation
    quality="hd"
)
```

## Example 3: Logo Design

### User Request
"Design a logo for a coffee roasting company called 'Mountain Peak Coffee'"

### Analysis
- **Type:** Logo design
- **Style:** Minimalist, professional
- **Format:** Square (versatile)
- **Model choice:** DALL-E 2 (multiple variations) → refine with DALL-E 3

### Prompt Strategy
Generate multiple concepts first, then refine the chosen one.

### Phase 1: Multiple Concepts (DALL-E 2)
```python
prompt = """Minimalist coffee company logo for 'Mountain Peak Coffee', geometric
mountain peak silhouette combined with coffee bean or cup element, earth tones
with deep brown and forest green, clean modern design, professional and organic
aesthetic, vector style, simple and memorable"""

response = client.images.generate(
    model="dall-e-2",
    prompt=prompt,
    size="1024x1024",
    n=4  # Generate 4 concepts
)

for idx, img in enumerate(response.data):
    img_data = requests.get(img.url).content
    Path(f"logo_concept_{idx+1}.png").write_bytes(img_data)
    print(f"Generated: logo_concept_{idx+1}.png")
```

### Phase 2: Refine Selected Concept (DALL-E 3)
```python
# After user selects concept #2
refined_prompt = """Professional coffee roasting company logo for 'Mountain Peak
Coffee', minimalist design with stylized mountain peak forming the outline of a
coffee cup, earth tone color palette with warm brown and sage green, clean
modern aesthetic, vector style, simple geometric shapes, suitable for packaging
and signage"""

response = client.images.generate(
    model="dall-e-3",
    prompt=refined_prompt,
    size="1024x1024",
    quality="hd"
)

final_logo = requests.get(response.data[0].url).content
Path("mountain_peak_coffee_logo_final.png").write_bytes(final_logo)
```

## Example 4: Product Photography

### User Request
"Generate product photos of wireless earbuds for an e-commerce listing"

### Analysis
- **Type:** Product photography
- **Style:** Clean, commercial
- **Use case:** E-commerce (need multiple angles)
- **Model choice:** Gemini 2.5 Flash (speed for multiple shots) or DALL-E 3

### Multiple Angle Generation

**Angle 1: Hero Shot**
```python
prompt1 = """Wireless earbuds with charging case, professional product
photography, white seamless background, studio lighting with soft shadows,
3/4 angle view showing earbuds next to open charging case, matte black finish
with subtle branding, clean minimal composition, high resolution, commercial
e-commerce quality"""

# Generate hero shot
```

**Angle 2: Detail Shot**
```python
prompt2 = """Wireless earbuds close-up detail shot, professional product
photography, white background, single earbud in focus showing ergonomic design
and touch controls, dramatic side lighting highlighting curves and texture,
commercial quality, high detail"""

# Generate detail shot
```

**Angle 3: In-Use Context**
```python
prompt3 = """Wireless earbuds lifestyle shot, person's hand holding single
earbud about to place in ear, natural lifestyle photography, clean white
background, professional lighting, modern and elegant, product focus with
shallow depth of field"""

# Generate lifestyle shot
```

### Implementation
```python
prompts = [
    "Wireless earbuds with charging case, professional product photography...",
    "Wireless earbuds close-up detail shot...",
    "Wireless earbuds lifestyle shot..."
]

for idx, prompt in enumerate(prompts, 1):
    response = client.models.generate_content(
        model="gemini-2.5-flash-image",
        contents=prompt,
        config=types.GenerateContentConfig(
            response_modalities=["TEXT", "IMAGE"],
            image_config=types.ImageConfig(aspect_ratio="1:1")
        )
    )

    for part in response.parts:
        if part.inline_data is not None:
            image = part.as_image()
            image.save(f"earbuds_angle_{idx}.png")
            print(f"Generated: earbuds_angle_{idx}.png")
```

## Example 5: Illustration for Children's Book

### User Request
"Create an illustration of a friendly dragon for a children's book"

### Analysis
- **Type:** Illustration
- **Style:** Friendly, colorful, age-appropriate
- **Use case:** Children's book (portrait orientation)
- **Model choice:** DALL-E 3 (excellent for illustration style)

### Prompt
```
Friendly cartoon dragon character for children's book, cute and approachable
design with big expressive eyes, rounded soft shapes, sitting pose with wings
folded, bright cheerful colors with purple and teal scales, children's book
illustration style, watercolor and digital art technique, whimsical and playful,
warm and inviting
```

### Implementation (DALL-E 3)
```python
prompt = """Friendly cartoon dragon character for children's book, cute with
big expressive eyes, rounded soft shapes, sitting pose, bright purple and teal
colors, children's book illustration style, watercolor technique, whimsical
and playful"""

response = client.images.generate(
    model="dall-e-3",
    prompt=prompt,
    size="1024x1792",  # Portrait for book page
    quality="hd"
)

image_data = requests.get(response.data[0].url).content
Path("friendly_dragon_illustration.png").write_bytes(image_data)

# DALL-E 3 provides the revised prompt it used
print(f"DALL-E's interpretation: {response.data[0].revised_prompt}")
```

## Example 6: Social Media Graphics

### User Request
"Create an Instagram post announcing a flash sale"

### Analysis
- **Type:** Social media graphic
- **Style:** Eye-catching, modern
- **Format:** Square (1:1 for Instagram)
- **Model choice:** Gemini 2.5 Flash (fast) or DALL-E 3

### Prompt
```
Instagram marketing graphic for flash sale announcement, bold modern design,
"FLASH SALE 50% OFF" text prominently displayed, vibrant gradient background
from coral to purple, shopping-related iconography (lightning bolt, shopping
bags), clean professional layout, energetic and exciting mood, social media
graphic design style
```

### Implementation (Gemini Pro - better with text)
```python
prompt = """Instagram marketing graphic for flash sale announcement, bold modern
design with "FLASH SALE 50% OFF" text prominently displayed in large bold font,
vibrant gradient background from coral pink to purple, lightning bolt and
shopping bag icons, clean professional layout, energetic and exciting mood,
social media graphic design style, text must be clear and legible"""

response = client.models.generate_content(
    model="gemini-3-pro-image-preview",
    contents=prompt,
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        image_config=types.ImageConfig(
            aspect_ratio="1:1",
            image_size="2K"
        )
    )
)

for part in response.parts:
    if part.inline_data is not None:
        image = part.as_image()
        image.save("flash_sale_instagram.png")
```

## Example 7: Concept Art

### User Request
"Generate concept art for a sci-fi video game character"

### Analysis
- **Type:** Concept art / character design
- **Style:** Detailed, professional concept art
- **Use case:** Game development reference
- **Model choice:** Gemini Pro (high detail) or DALL-E 3 (creative)

### Prompt
```
Sci-fi video game character concept art, female space marine in advanced combat
armor, tactical military aesthetic with futuristic tech elements, armored
exosuit with glowing blue energy panels, helmet under arm showing determined
expression, dynamic standing pose, detailed digital concept art style, gray and
blue color scheme with orange accent lights, professional game art quality,
character design sheet format with front view
```

### Implementation (DALL-E 3)
```python
prompt = """Sci-fi video game character concept art, female space marine in
advanced combat armor, futuristic tactical exosuit with glowing blue energy
panels, determined expression, dynamic pose, detailed digital concept art,
gray and blue with orange accents, professional game art quality"""

response = client.images.generate(
    model="dall-e-3",
    prompt=prompt,
    size="1024x1792",  # Portrait for character
    quality="hd"
)

concept_art = requests.get(response.data[0].url).content
Path("scifi_character_concept.png").write_bytes(concept_art)
```

## Example 8: Interior Design Visualization

### User Request
"Show me what a modern minimalist living room would look like"

### Analysis
- **Type:** Interior design visualization
- **Style:** Photorealistic, modern
- **Use case:** Design inspiration
- **Model choice:** Gemini Pro (high detail) or DALL-E 3

### Prompt
```
Modern minimalist living room interior, spacious open layout with floor-to-ceiling
windows, natural light flooding the space, neutral color palette with white walls
and light wood floors, low-profile gray sofa, simple coffee table, indoor plants,
minimal decor, clean lines, Scandinavian design influence, architectural interior
photography, bright and airy atmosphere, high detail
```

### Implementation (Gemini)
```python
prompt = """Modern minimalist living room interior, spacious with floor-to-ceiling
windows, natural light, neutral colors with white walls and light wood floors,
gray sofa, simple furniture, indoor plants, Scandinavian design, architectural
photography, bright and airy, high detail"""

response = client.models.generate_content(
    model="gemini-3-pro-image-preview",
    contents=prompt,
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"],
        image_config=types.ImageConfig(
            aspect_ratio="16:9",
            image_size="4K"
        )
    )
)
```

## Example 9: Abstract Art

### User Request
"Create abstract art representing the concept of innovation"

### Analysis
- **Type:** Abstract art
- **Style:** Artistic, conceptual
- **Use case:** Office decor / presentation
- **Model choice:** DALL-E 3 (excellent for artistic interpretation)

### Prompt
```
Abstract art representing innovation and creativity, dynamic flowing shapes and
geometric forms interconnecting, vibrant color palette with electric blue,
orange, and purple, sense of movement and energy, digital abstract art style,
modern and inspiring, suitable for contemporary office space, explosive and
forward-thinking composition
```

### Implementation (DALL-E 3)
```python
prompt = """Abstract art representing innovation and creativity, dynamic flowing
shapes and geometric forms, vibrant electric blue, orange and purple colors,
sense of movement and energy, digital abstract art, modern and inspiring"""

response = client.images.generate(
    model="dall-e-3",
    prompt=prompt,
    size="1024x1024",
    quality="hd"
)
```

## Example 10: Food Photography

### User Request
"Generate an appetizing photo of a gourmet burger for our restaurant menu"

### Analysis
- **Type:** Food photography
- **Style:** Appetizing, professional
- **Use case:** Menu/marketing
- **Model choice:** Gemini Pro or DALL-E 3 (high quality needed)

### Prompt
```
Gourmet burger professional food photography, artisan sesame bun with perfectly
grilled beef patty, melted aged cheddar, crisp lettuce, heirloom tomato slices,
caramelized onions, on rustic wooden board, side of golden french fries, natural
lighting with soft shadows, 45-degree overhead angle, appetizing and fresh,
restaurant quality food photography, steam rising, high detail
```

### Implementation (DALL-E 3)
```python
prompt = """Gourmet burger professional food photography, artisan bun with
grilled beef patty, melted cheddar, fresh vegetables, rustic wooden board,
golden fries on side, natural lighting, overhead angle, appetizing and fresh,
restaurant quality, high detail"""

response = client.images.generate(
    model="dall-e-3",
    prompt=prompt,
    size="1024x1024",
    quality="hd"
)

burger_photo = requests.get(response.data[0].url).content
Path("gourmet_burger_menu.png").write_bytes(burger_photo)
```

## Best Practices from Examples

### 1. Always Include:
- Subject/main focus
- Style (photography, illustration, art style)
- Lighting and mood
- Technical details (lens, composition)
- Use case context

### 2. Model Selection:
- **Speed needed** → Gemini Flash or DALL-E 2
- **Quality needed** → Gemini Pro or DALL-E 3
- **Multiple variations** → DALL-E 2
- **Text in image** → Gemini Pro

### 3. Format Selection:
- **Square (1:1)** → Social media, logos, portraits
- **Portrait (3:4, 9:16)** → Characters, mobile, stories
- **Landscape (16:9, 21:9)** → Headers, scenery, wide shots

### 4. Iteration:
- Start with clear base prompt
- Generate and evaluate
- Refine specific elements
- Use chat/conversation for Gemini refinements

### 5. Saving and Organization:
- Use descriptive filenames
- Save in appropriate resolution
- Keep track of prompts used
- Organize by project/category

## Common Patterns Summary

**Photography:** `[Subject] + professional [type] photography + [lighting] + [camera specs] + high detail`

**Illustration:** `[Subject] + [art style] illustration + [color palette] + [mood] + [technique]`

**Design:** `[Concept] + [style] design + [elements] + [colors] + professional + [format]`

**Art:** `[Concept/subject] + [art movement/style] + [mood/feeling] + [technical approach]`

Use these examples as templates and adapt to your specific needs!
